/*
 * 
 * “Commons Clause” License Condition v1.0
 * 
 * The Software is provided to you by the Licensor under the License, 
 * as defined below, subject to the following condition.
 * 
 * Without limiting other conditions in the License, the grant of rights under the License 
 * will not include, and the License does not grant to you, the right to Sell the Software.
 * 
 * For purposes of the foregoing, “Sell” means practicing any or all of the rights granted 
 * to you under the License to provide to third parties, for a fee or other consideration 
 * (including without limitation fees for hosting or consulting/ support services related to 
 * the Software), a product or service whose value derives, entirely or substantially, from the 
 * functionality of the Software. Any license notice or attribution required by the License 
 * must also include this Commons Clause License Condition notice.
 * 
 * Software: genestrip
 * 
 * License: Apache 2.0
 * 
 * Licensor: Daniel Pfeifer (daniel.pfeifer@progotec.de)
 * 
 */
package org.metagene.genestrip.match;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;

import org.metagene.genestrip.fastq.AbstractFastqReader;
import org.metagene.genestrip.match.ResultReporter.Result;
import org.metagene.genestrip.match.ResultReporter.StatsPerTaxid;
import org.metagene.genestrip.store.KMerSortedArray;
import org.metagene.genestrip.store.KMerUniqueCounter;
import org.metagene.genestrip.util.ByteArrayUtil;
import org.metagene.genestrip.util.CGAT;
import org.metagene.genestrip.util.StreamProvider;
import org.metagene.genestrip.util.StreamProvider.ByteCountingInputStreamAccess;

import it.unimi.dsi.fastutil.objects.Object2IntMap;

public class FastqKMerMatcher extends AbstractFastqReader {
	public static final long DEFAULT_LOG_UPDATE_CYCLE = 1000000;

	private final KMerSortedArray<String> kmerStore;

	protected KMerUniqueCounter uniqueCounter;
	protected TaxidStatsTrie root;

	private ByteCountingInputStreamAccess byteCountAccess;
	private long fastqFileSize;
	private long startTime;
	private long indexedC;

	private OutputStream indexed;

	private long logUpdateCycle = DEFAULT_LOG_UPDATE_CYCLE;

	// A PrintStream is implicitly synchronized. So we don't need to worry about
	// multi threading when using it.
	protected PrintStream out;

	public FastqKMerMatcher(KMerSortedArray<String> kmerStore, int maxReadSize, int maxQueueSize, int consumerNumber) {
		super(kmerStore.getK(), maxReadSize, maxQueueSize, consumerNumber);
		this.kmerStore = kmerStore;
	}

	@Override
	protected ReadEntry createReadEntry(int maxReadSizeBytes) {
		return new MyReadEntry(maxReadSizeBytes, k);
	}

	public Result runClassifier(File fastq, File filteredFile, File krakenOutStyleFile,
			KMerUniqueCounter uniqueCounter) throws IOException {
		byteCountAccess = StreamProvider.getByteCountingInputStreamForFile(fastq, false);
		fastqFileSize = Files.size(fastq.toPath());

		if (filteredFile != null) {
			indexed = StreamProvider.getOutputStreamForFile(filteredFile);
		} else {
			indexed = null;
		}
		if (krakenOutStyleFile != null) {
			// A PrintStream is implicitly synchronized. So we don't need to worry about
			// multi threading when using it.
			out = new PrintStream(StreamProvider.getOutputStreamForFile(krakenOutStyleFile));
		} else {
			out = null;
		}

		startTime = System.currentTimeMillis();

		initRoot();
		this.uniqueCounter = uniqueCounter;
		if (uniqueCounter != null) {
			uniqueCounter.clear();
		}

		readFastq(byteCountAccess.getInputStream());

		if (indexed != null) {
			indexed.close();
			out.close();
		}
		byteCountAccess.getInputStream().close();

		List<StatsPerTaxid> allStats = new ArrayList<StatsPerTaxid>();
		root.collectValues(allStats);
		root = null;

		if (uniqueCounter != null) {
			Object2IntMap<String> counts = uniqueCounter.getUniqueKmerCounts();
			for (StatsPerTaxid stats : allStats) {
				stats.uniqueKmers = counts.getInt(stats.getTaxid());
			}
			this.uniqueCounter = null;
		} else {
			for (StatsPerTaxid stats : allStats) {
				stats.uniqueKmers = -1;
			}
		}

		return new Result(allStats, reads, kMers);
	}

	protected void initRoot() {
		root = new TaxidStatsTrie();
	}

	@Override
	protected void nextEntry(ReadEntry entry) throws IOException {
		MyReadEntry myEntry = (MyReadEntry) entry;
		myEntry.bufferPos = 0;
		myEntry.readTaxId = null;

		boolean found = classifyRead(myEntry, false);
		if (!found) {
			found = classifyRead(myEntry, true);
		}
		if (found) {
			if (indexed != null) {
				rewriteInput(entry, indexed);
			}
			if (out != null) {
				synchronized (out) {
					myEntry.flush(out);
				}
			}
		}
	}

	@Override
	protected void start() throws IOException {
		indexedC = 0;
	}

	@Override
	protected void updateWriteStats() {
		indexedC++;
	}

	@Override
	protected void log() {
		if (reads % logUpdateCycle == 0) {
			if (logger.isInfoEnabled()) {
				double ratio = byteCountAccess.getBytesRead() / (double) fastqFileSize;
				long stopTime = System.currentTimeMillis();

				double diff = (stopTime - startTime);
				double totalTime = diff / ratio;
				double totalHours = totalTime / 1000 / 60 / 60;

				logger.info("Elapsed hours:" + diff / 1000 / 60 / 60);
				logger.info("Estimated total hours:" + totalHours);
				logger.info("Reads processed: " + reads);
				logger.info("Indexed: " + indexedC);
				logger.info("Indexed ratio:" + ((double) indexedC) / reads);
			}
		}
	}

	public long getLogUpdateCycle() {
		return logUpdateCycle;
	}

	public void setLogUpdateCycle(long logUpdateCycle) {
		this.logUpdateCycle = logUpdateCycle;
	}

	@Override
	protected void done() throws IOException {
		if (logger.isInfoEnabled()) {
			logger.info("Total indexed: " + indexedC);
			long stopTime = System.currentTimeMillis();
			double diff = (stopTime - startTime);
			logger.info("Elapsed total ms:" + diff);
			logger.info("Read file time ms: " + getMillis());
		}
	}

	protected boolean classifyRead(MyReadEntry entry, boolean reverse) {
		boolean found = false;
		int prints = 0;
		int taxIdCounter = 0;

		String taxid = null;
		int max = entry.readSize - k;
		String lastTaxid = null;
		int contigLen = 0;
		StatsPerTaxid stats = null;

		for (int i = 0; i <= max; i++) {
			long kmer = reverse ? CGAT.kMerToLongReverse(entry.read, i, k, null)
					: CGAT.kMerToLongStraight(entry.read, i, k, null);

			taxid = kmerStore.getLong(kmer, entry.indexPos);
			if (contigLen > 0 && taxid != lastTaxid) {
				if (out != null) {
					printKrakenStyleOut(entry, lastTaxid, contigLen, prints++, reverse);
				}
				if (stats != null) {
					synchronized (stats) {
						stats.contigs++;
						if (contigLen > stats.maxContigLen) {
							stats.maxContigLen = contigLen;
						}
					}
				}
				contigLen = 0;
			}
			contigLen++;
			lastTaxid = taxid;
			if (taxid != null) {
				stats = root.get(taxid);
				if (stats == null) {
					stats = root.create(taxid);
				}
				synchronized (stats) {
					stats.kmers++;
					// TODO: This is wrong if there are two taxids (or more) in the same read. The
					// read counter is only increased for the first one.
					if (!found) {
						stats.reads++;
					}
					if (entry.readTaxId != taxid) {
						taxIdCounter++;
					}
					if (taxIdCounter == 1) {
						entry.readTaxId = taxid;
					}
					found = true;
				}
				if (uniqueCounter != null) {
					uniqueCounter.put(kmer, taxid, entry.indexPos[0]);
				}
			} else {
				stats = null;
			}
		}
		if (found) {
			if (contigLen > 0) {
				if (out != null) {
					printKrakenStyleOut(entry, lastTaxid, contigLen, prints, reverse);
				}
				if (stats != null) {
					synchronized (stats) {
						stats.contigs++;
						if (contigLen > stats.maxContigLen) {
							stats.maxContigLen = contigLen;
						}
					}
				}
			}
			if (taxIdCounter > 1) {
				entry.readTaxId = null;
			}
		}

		return found;
	}

	protected void printKrakenStyleOut(MyReadEntry entry, String taxid, int contigLen, int state, boolean reverse) {
		if (state != 0) {
			entry.printChar(' ');
		}
		if (taxid == null) {
			entry.printChar('0');
		} else {
			entry.printString(taxid);
		}
		entry.printChar(':');
		entry.printInt(contigLen);
	}

	public static class MyReadEntry extends ReadEntry {
		public final byte[] buffer;
		public int bufferPos;
		public String readTaxId;
		public long[] indexPos;

		public MyReadEntry(int maxReadSizeBytes, int k) {
			super(maxReadSizeBytes);

			buffer = new byte[maxReadSizeBytes];
			indexPos = new long[1];
		}

		public void printChar(char c) {
			buffer[bufferPos++] = (byte) c;
		}

		public void printString(String s) {
			int len = s.length();
			for (int i = 0; i < len; i++) {
				buffer[bufferPos++] = (byte) s.charAt(i);
			}
		}

		public void printBytes(byte[] bytes) {
			for (int i = 0; i < bytes.length; i++) {
				byte b = bytes[i];
				if (b == 0) {
					return;
				}
				buffer[bufferPos++] = b;
			}
		}

		public void printInt(int value) {
			bufferPos = ByteArrayUtil.intToByteArray(value, buffer, bufferPos);
		}

		public void flush(PrintStream out) {
			if (readTaxId != null) {
				out.print("U\t");
			} else {
				out.print("C\t");
			}
			ByteArrayUtil.print(readDescriptor, out);
			out.print('\t');
			if (readTaxId != null) {
				out.print(readTaxId);
			} else {
				out.print('0');
			}
			out.print('\t');
			out.print(readSize);
			out.print('\t');
			out.write(buffer, 0, bufferPos);
			out.println();
		}
	}

	protected static class TaxidStatsTrie {
		private TaxidStatsTrie[] children;
		private StatsPerTaxid value;

		public StatsPerTaxid get(String taxid) {
			int index;
			int end = taxid.length();
			TaxidStatsTrie node = this, child;
			for (int i = 0; i < end; i++, node = child) {
				index = taxid.charAt(i) - '0';
				if (index > 9 || index < 0) {
					return null;
				}
				if (node.children == null) {
					return null;
				}
				child = node.children[index];
				if (child == null) {
					return null;
				}
			}
			return node.value;
		}

		public synchronized StatsPerTaxid create(String taxid) {
			int index;
			int end = taxid.length();
			TaxidStatsTrie node = this, child;
			for (int i = 0; i < end; i++, node = child) {
				index = taxid.charAt(i) - '0';
				if (index > 9 || index < 0) {
					return null;
				}
				if (node.children == null) {
					node.children = new TaxidStatsTrie[10];
				}
				child = node.children[index];
				if (child == null) {
					child = new TaxidStatsTrie();
					node.children[index] = child;
				}
			}
			if (node.value == null) {
				node.value = new StatsPerTaxid(taxid);
			}
			return node.value;
		}

		public void collectValues(List<StatsPerTaxid> list) {
			if (value != null) {
				list.add(value);
			}
			if (children != null) {
				for (int i = 0; i < 10; i++) {
					if (children[i] != null) {
						children[i].collectValues(list);
					}
				}
			}
		}
	}

}
