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
package org.metagene.genestrip.store;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;

import org.metagene.genestrip.fastq.AbstractFastqReader;
import org.metagene.genestrip.util.ByteArrayUtil;
import org.metagene.genestrip.util.CGAT;
import org.metagene.genestrip.util.StreamProvider;
import org.metagene.genestrip.util.StreamProvider.ByteCountingInputStreamAccess;

import it.unimi.dsi.fastutil.longs.Long2ObjectMap;
import it.unimi.dsi.fastutil.longs.Long2ObjectOpenHashMap;
import it.unimi.dsi.fastutil.objects.ObjectCollection;

public class FastqKMerMatcher extends AbstractFastqReader {
	public static final long DEFAULT_LOG_UPDATE_CYCLE = 1000000;

	private final KMerStore<String> kmerStore;
	protected final KmerDuplicationCount duplicationCount;

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

	public FastqKMerMatcher(KMerStore<String> kmerStore, int maxReadSize, int maxQueueSize, int consumerNumber,
			boolean withDupCount) {
		super(kmerStore.getK(), maxReadSize, maxQueueSize, consumerNumber);
		this.kmerStore = kmerStore;
		duplicationCount = withDupCount ? new KmerDuplicationCount(k) : null;
	}

	@Override
	protected ReadEntry createReadEntry(int maxReadSizeBytes) {
		return new MyReadEntry(maxReadSizeBytes, k);
	}

	public Result runClassifier(File fastq, File filteredFile, File krakenOutStyleFile)
			throws IOException {
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
		if (duplicationCount != null) {
			duplicationCount.clear();
		}

		readFastq(byteCountAccess.getInputStream());

		if (indexed != null) {
			indexed.close();
			out.close();
		}
		byteCountAccess.getInputStream().close();

		List<StatsPerTaxid> counts = new ArrayList<FastqKMerMatcher.StatsPerTaxid>();
		root.collectValues(counts);

		if (duplicationCount != null) {
			for (StatsPerTaxid stats : counts) {
				stats.uniqueKmers = duplicationCount.getUniqeKmers(stats.taxid);
				if (stats.kmers != duplicationCount.getKMerCount(stats.taxid)) {
					if (logger.isWarnEnabled()) {
						logger.warn("Inconsistent kmer counts for taxid: " + stats.taxid);
					}
				}
			}
		} else {
			for (StatsPerTaxid stats : counts) {
				stats.uniqueKmers = -1;
			}
		}
		
		return new Result(counts, reads, kMers);
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
			taxid = kmerStore.get(entry.read, i, reverse);
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
				if (duplicationCount != null) {
					if (!duplicationCount.put(taxid, entry, i, reverse)) {
						if (logger.isInfoEnabled()) {
							synchronized (logger) {
								logger.warn("Updating duplication count failed for kmer under taxid " + taxid);
							}
						}
					}
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

		public MyReadEntry(int maxReadSizeBytes, int k) {
			super(maxReadSizeBytes);

			buffer = new byte[maxReadSizeBytes];
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

	public static class KmerDuplicationCount {
		private final Long2ObjectMap<Entry> map;
		private final int k;

		public KmerDuplicationCount(int k) {
			this.k = k;
			map = new Long2ObjectOpenHashMap<Entry>();
		}

		public void clear() {
			map.clear();
		}

		protected boolean put(String taxid, MyReadEntry entry, int start, boolean reverse) {
			long kmer;
			if (reverse) {
				kmer = CGAT.kMerToLongReverse(entry.read, start, k, null);
			} else {
				kmer = CGAT.kMerToLongStraight(entry.read, start, k, null);
			}

			synchronized (map) {
				Entry e = map.get(kmer);
				if (e == null) {
					e = new Entry();
					e.taxid = taxid;
					map.put(kmer, e);
				} else if (e.taxid != taxid) {
					return false;
				}
				e.count++;
			}

			return true;
		}

		public int getUniqeKmers(String taxid) {
			int res = 0;
			ObjectCollection<Entry> entries = map.values();
			for (Entry e : entries) {
				if (e.taxid == taxid) {
					res++;
				}
			}
			return res;
		}

		public int getKMerCount(String taxid) {
			int res = 0;
			ObjectCollection<Entry> entries = map.values();
			for (Entry e : entries) {
				if (e.taxid == taxid) {
					res += e.count;
				}
			}
			return res;
		}

		public int getKmerCount(long kmerCode) {
			Entry e = map.get(kmerCode);
			return e == null ? 0 : e.count;
		}

		private static class Entry {
			public int count;
			public String taxid;
		}
	}

	public static class StatsPerTaxid {
		protected String taxid;
		protected long reads;
		protected long uniqueKmers;
		protected long kmers;
		protected int maxContigLen;
		protected int contigs;

		public StatsPerTaxid(String taxid) {
			this.taxid = taxid;
		}

		public String getTaxid() {
			return taxid;
		}

		public int getContigs() {
			return contigs;
		}

		public long getKMers() {
			return kmers;
		}

		public int getMaxContigLen() {
			return maxContigLen;
		}

		public long getReads() {
			return reads;
		}

		public long getUniqueKMers() {
			return uniqueKmers;
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
	
	public static class Result {
		private final List<StatsPerTaxid> stats;
		private final long totalReads;
		private final long totalKMers;
		
		public Result(List<StatsPerTaxid> stats, long totalReads, long totalKMers) {
			this.stats = stats;
			this.totalReads = totalReads;
			this.totalKMers = totalKMers;
		}
		
		public List<StatsPerTaxid> getStats() {
			return stats;
		}
		
		public long getTotalKMers() {
			return totalKMers;
		}
		
		public long getTotalReads() {
			return totalReads;
		}
	}
}