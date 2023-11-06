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
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.metagene.genestrip.fastq.AbstractFastqReader;
import org.metagene.genestrip.io.StreamProvider;
import org.metagene.genestrip.io.StreamProvider.ByteCountingInputStreamAccess;
import org.metagene.genestrip.store.KMerSortedArray;
import org.metagene.genestrip.store.KMerUniqueCounter;
import org.metagene.genestrip.store.KMerUniqueCounterBits;
import org.metagene.genestrip.tax.TaxTree;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;
import org.metagene.genestrip.util.ByteArrayUtil;
import org.metagene.genestrip.util.CGAT;
import org.metagene.genestrip.util.DigitTrie;

import it.unimi.dsi.fastutil.objects.Object2LongMap;

public class FastqKMerMatcher2 extends AbstractFastqReader {
	public static final long DEFAULT_LOG_UPDATE_CYCLE = 1000000;

	private final KMerSortedArray<TaxIdNode> kmerStore;
	private final int maxKmerResCounts;

	protected KMerUniqueCounter uniqueCounter;
	protected TaxidStatsTrie root;

	private ByteCountingInputStreamAccess byteCountAccess;
	private final int paths;
	private int coveredCounter;
	private int totalCount;
	private File currentFastq;
	private long fastqsFileSize;
	private long coveredFilesSize;
	private long startTime;
	private long indexedC;
	private final TaxTree taxTree;
	protected final double maxReadTaxErrorCount;

	private OutputStream indexed;

	private final Integer maxReadSize;

	private long logUpdateCycle = DEFAULT_LOG_UPDATE_CYCLE;

	// A PrintStream is implicitly synchronized. So we don't need to worry about
	// multi-threading when using it.
	protected PrintStream out;

	public FastqKMerMatcher2(KMerSortedArray<TaxIdNode> kmerStore, int maxReadSize, int maxQueueSize,
			int consumerNumber, int maxKmerResCounts, TaxTree taxTree, double maxReadTaxErrorCount) {
		super(kmerStore.getK(), maxReadSize, maxQueueSize, consumerNumber);
		this.kmerStore = kmerStore;
		this.maxReadSize = maxReadSize;
		this.maxKmerResCounts = maxKmerResCounts;
		this.taxTree = taxTree;
		this.maxReadTaxErrorCount = maxReadTaxErrorCount;
		paths = 10; // TODO: Make configurable
		taxTree.initCountSize(consumerNumber == 0 ? 1 : consumerNumber);
	}

	@Override
	protected ReadEntry createReadEntry(int maxReadSizeBytes) {
		return new MyReadEntry(maxReadSizeBytes, out != null, paths);
	}

	public MatchingResult runMatcher(File fastq, File filteredFile, File krakenOutStyleFile,
			KMerUniqueCounter uniqueCounter) throws IOException {
		return runMatcher(Collections.singletonList(fastq), filteredFile, krakenOutStyleFile, uniqueCounter);
	}

	public MatchingResult runMatcher(List<File> fastqs, File filteredFile, File krakenOutStyleFile,
			KMerUniqueCounter uniqueCounter) throws IOException {
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
		initUniqueCounter(uniqueCounter);

		totalCount = fastqs.size();
		for (File fastq : fastqs) {
			fastqsFileSize += Files.size(fastq.toPath());
		}
		coveredCounter = 0;
		for (File fastq : fastqs) {
			currentFastq = fastq;
			byteCountAccess = StreamProvider.getByteCountingInputStreamForFile(fastq, false);
			readFastq(byteCountAccess.getInputStream());
			coveredFilesSize += byteCountAccess.getBytesRead();
			byteCountAccess.getInputStream().close();
			coveredCounter++;
		}

		if (indexed != null) {
			indexed.close();
			out.close();
		}

		List<CountsPerTaxid> allStats = new ArrayList<CountsPerTaxid>();
		root.collect(allStats);
		root = null;
		Map<String, CountsPerTaxid> taxid2Stats = new HashMap<String, CountsPerTaxid>();
		for (CountsPerTaxid stats : allStats) {
			taxid2Stats.put(stats.getTaxid(), stats);
		}

		Map<String, short[]> countMap = null;
		if (uniqueCounter != null) {
			Object2LongMap<String> counts = uniqueCounter.getUniqueKmerCounts();
			for (CountsPerTaxid stats : allStats) {
				stats.uniqueKmers = counts.getLong(stats.getTaxid());
			}
			if (uniqueCounter instanceof KMerUniqueCounterBits) {
				if (((KMerUniqueCounterBits) uniqueCounter).isWithCounts()) {
					countMap = ((KMerUniqueCounterBits) uniqueCounter).getMaxCountsCounts(maxKmerResCounts);
					for (CountsPerTaxid stats : allStats) {
						stats.maxKMerCounts = countMap.get(stats.getTaxid());
					}
				}
			}
			this.uniqueCounter = null;
		} else {
			for (CountsPerTaxid stats : allStats) {
				stats.uniqueKmers = -1;
			}
		}

		return new MatchingResult(kmerStore.getK(), taxid2Stats, reads, kMers,
				countMap == null ? null : countMap.get(null));
	}

	protected void initRoot() {
		root = new TaxidStatsTrie();
	}

	protected void initUniqueCounter(KMerUniqueCounter uniqueCounter) {
		this.uniqueCounter = uniqueCounter;
		if (uniqueCounter != null) {
			uniqueCounter.clear();
		}
	}

	@Override
	protected void nextEntry(ReadEntry entry, int index) throws IOException {
		MyReadEntry myEntry = (MyReadEntry) entry;
		myEntry.bufferPos = 0;
		myEntry.readTaxErrorCount = 0;

		myEntry.usedPaths = 0;
		myEntry.classNode = null;
		for (int i = 0; i < paths; i++) {
			myEntry.readTaxIdNode[i] = null;
			myEntry.counts[i] = 0;
		}

		boolean found = matchRead(myEntry, index, false);
		if (!found) {
			found = matchRead(myEntry, index, true);
		}
		afterMatch(myEntry, found);
	}

	protected void afterMatch(MyReadEntry myEntry, boolean found) throws IOException {
		if (found) {
			if (indexed != null) {
				rewriteInput(myEntry, indexed);
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
		if (logger.isInfoEnabled()) {
			logger.info("Processing fastq file (" + (coveredCounter + 1) + "/" + totalCount + "): " + currentFastq);
		}
	}

	@Override
	protected void updateWriteStats() {
		indexedC++;
	}

	@Override
	protected void log() {
		if (reads % logUpdateCycle == 0) {
			if (logger.isInfoEnabled()) {
				double ratio = (coveredFilesSize + byteCountAccess.getBytesRead()) / (double) fastqsFileSize;
				long stopTime = System.currentTimeMillis();

				double diff = (stopTime - startTime);
				double totalTime = diff / ratio;
				double totalHours = totalTime / 1000 / 60 / 60;

				logger.info("Elapsed hours: " + diff / 1000 / 60 / 60);
				logger.info("Estimated total hours: " + totalHours);
				logger.info("Reads processed: " + reads);
				logger.info("Indexed: " + indexedC);
				logger.info("Indexed ratio: " + ((double) indexedC) / reads);
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
			logger.info("Elapsed total ms: " + diff);
			logger.info("Read file time ms: " + getMillis());
		}
	}

	protected boolean matchRead(MyReadEntry entry, int index, boolean reverse) {
		boolean found = false;
		int prints = 0;

		TaxIdNode taxid = null;
		int max = entry.readSize - k;
		TaxIdNode lastTaxid = null;
		int contigLen = 0;
		CountsPerTaxid stats = null;

		for (int i = 0; i <= max; i++) {
			long kmer = reverse ? CGAT.kMerToLongReverse(entry.read, i, k, null)
					: CGAT.kMerToLongStraight(entry.read, i, k, null);

			taxid = kmer == -1 ? null : kmerStore.getLong(kmer, entry.indexPos);
			if (taxid == null) {
				entry.readTaxErrorCount++;
				if ((maxReadTaxErrorCount >= 1 && entry.readTaxErrorCount <= maxReadTaxErrorCount)
						|| (entry.readTaxErrorCount <= maxReadTaxErrorCount * max)) {
					
				}
			} else {
				updateReadTaxid(taxid, entry, index);
			}
			if (taxid != lastTaxid) {
				if (contigLen > 0) {
					if (out != null) {
						printKrakenStyleOut(entry, lastTaxid, contigLen, prints++, reverse);
					}
					if (stats != null) {
						synchronized (stats) {
							stats.contigs++;
							if (contigLen > stats.maxContigLen) {
								stats.maxContigLen = contigLen;
								for (int j = 0; j < entry.readSize; j++) {
									stats.maxContigDescriptor[j] = entry.readDescriptor[j];
								}
							}
						}
					}
					contigLen = 0;
				}
			}
			contigLen++;
			lastTaxid = taxid;
			if (taxid != null) {
				stats = root.get(taxid.getTaxId());
				if (stats == null) {
					synchronized (root) {
						stats = root.get(taxid.getTaxId(), maxReadSize);
					}
				}
				synchronized (stats) {
					stats.kmers++;
					found = true;
				}
				if (uniqueCounter != null) {
					uniqueCounter.put(kmer, taxid.getTaxId(), entry.indexPos[0]);
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
							for (int j = 0; j < entry.readSize; j++) {
								stats.maxContigDescriptor[j] = entry.readDescriptor[j];
							}
							stats.maxContigDescriptor[entry.readSize] = 0;
						}
					}
				}
			}
			int ties = 0;
			for (int i = 0; i < entry.usedPaths; i++) {
				short sum = taxTree.sumCounts(entry.readTaxIdNode[i], index, entry);
				if (sum >= entry.counts[0]) {
					entry.counts[0] = sum;
					entry.readTaxIdNode[0] = entry.readTaxIdNode[i];
					ties = 0;
				} else if (sum == entry.counts[0]) {
					ties++;
					entry.counts[ties] = sum;
					entry.readTaxIdNode[ties] = entry.readTaxIdNode[i];
				}
			}
			TaxIdNode node = entry.readTaxIdNode[0];
			for (int i = 1; i <= ties; i++) {
				node = taxTree.getLeastCommonAncestor(node, entry.readTaxIdNode[i]);
			}
			entry.classNode = node;
		}

		return found;
	}

	protected void updateReadTaxid(TaxIdNode node, MyReadEntry entry, int index) {
		if (node != null) {
			taxTree.incCount(node, index, entry);
		}
		boolean found = false;
		for (int i = 0; i < entry.usedPaths; i++) {
			if (taxTree.isAncestorOf(node, entry.readTaxIdNode[i])) {
				entry.readTaxIdNode[i] = node;
				found = true;
				break;
			} else if (taxTree.isAncestorOf(node, entry.readTaxIdNode[i])) {
				found = true;
				break;
			}
		}
		if (!found) {
			if (entry.usedPaths < paths) {
				entry.readTaxIdNode[entry.usedPaths] = node;
				entry.usedPaths++;
			}
		}
	}

	protected void printKrakenStyleOut(MyReadEntry entry, TaxIdNode taxid, int contigLen, int state, boolean reverse) {
		if (state != 0) {
			entry.printChar(' ');
		}
		if (taxid == null) {
			entry.printChar('0');
		} else {
			entry.printString(taxid.getTaxId());
		}
		entry.printChar(':');
		entry.printInt(contigLen);
	}

	public static class MyReadEntry extends ReadEntry {
		public final byte[] buffer;
		public int bufferPos;

		public int usedPaths;
//		public String readTaxId;
		public TaxIdNode[] readTaxIdNode;
		public short counts[];
		public int readTaxErrorCount;
		public long[] indexPos;
		public TaxIdNode classNode;

		public MyReadEntry(int maxReadSizeBytes, boolean enablePrint, int paths) {
			super(maxReadSizeBytes);

			buffer = enablePrint ? new byte[maxReadSizeBytes * 4] : null; // It has to be rather long in some cases...
			readTaxIdNode = new TaxIdNode[paths];
			counts = new short[paths];
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
			if (classNode == null) {
				out.print("U\t");
			} else {
				out.print("C\t");
			}
			ByteArrayUtil.print(readDescriptor, out);
			out.print('\t');
			if (classNode == null) {
				out.print('0');
			} else {
				out.print(classNode.getTaxId());
			}
			out.print('\t');
			out.print(readSize);
			out.print('\t');
			out.write(buffer, 0, bufferPos);
			out.println();
		}
	}

	protected static class TaxidStatsTrie extends DigitTrie<CountsPerTaxid> {
		@Override
		protected CountsPerTaxid createInGet(String digits, Object createContext) {
			return new CountsPerTaxid(digits, (int) createContext);
		}
	}
}
