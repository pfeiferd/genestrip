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
import java.io.Serializable;
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
import org.metagene.genestrip.util.ByteArrayUtil;
import org.metagene.genestrip.util.CGAT;
import org.metagene.genestrip.util.DigitTrie;

import it.unimi.dsi.fastutil.objects.Object2LongMap;

public class FastqKMerMatcher extends AbstractFastqReader {
	public static final long DEFAULT_LOG_UPDATE_CYCLE = 1000000;

	private final KMerSortedArray<String> kmerStore;

	protected KMerUniqueCounter uniqueCounter;
	protected TaxidStatsTrie root;

	private ByteCountingInputStreamAccess byteCountAccess;
	private int coveredCounter;
	private int totalCount;
	private File currentFastq;
	private long fastqsFileSize;
	private long coveredFilesSize;
	private long startTime;
	private long indexedC;

	private OutputStream indexed;
	
	private final Integer maxReadSize;

	private long logUpdateCycle = DEFAULT_LOG_UPDATE_CYCLE;

	// A PrintStream is implicitly synchronized. So we don't need to worry about
	// multi threading when using it.
	protected PrintStream out;

	public FastqKMerMatcher(KMerSortedArray<String> kmerStore, int maxReadSize, int maxQueueSize, int consumerNumber) {
		super(kmerStore.getK(), maxReadSize, maxQueueSize, consumerNumber);
		this.kmerStore = kmerStore;
		this.maxReadSize = maxReadSize;
	}

	@Override
	protected ReadEntry createReadEntry(int maxReadSizeBytes) {
		return new MyReadEntry(maxReadSizeBytes, out != null);
	}

	public Result runClassifier(File fastq, File filteredFile, File krakenOutStyleFile, KMerUniqueCounter uniqueCounter)
			throws IOException {
		return runClassifier(Collections.singletonList(fastq), filteredFile, krakenOutStyleFile, uniqueCounter);
	}

	public Result runClassifier(List<File> fastqs, File filteredFile, File krakenOutStyleFile,
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

		List<StatsPerTaxid> allStats = new ArrayList<StatsPerTaxid>();
		root.collect(allStats);
		root = null;
		Map<String, StatsPerTaxid> taxid2Stats = new HashMap<String, FastqKMerMatcher.StatsPerTaxid>();
		for (StatsPerTaxid stats : allStats) {
			taxid2Stats.put(stats.getTaxid(), stats);			
		}

		Map<String, short[]> countMap = null;
		if (uniqueCounter != null) {
			Object2LongMap<String> counts = uniqueCounter.getUniqueKmerCounts();
			for (StatsPerTaxid stats : allStats) {
				stats.uniqueKmers = counts.getLong(stats.getTaxid());
			}
			if (uniqueCounter instanceof KMerUniqueCounterBits) {
				if (((KMerUniqueCounterBits) uniqueCounter).isWithCounts()) {
					countMap = ((KMerUniqueCounterBits) uniqueCounter).getMaxCountsCounts(200);
					for (StatsPerTaxid stats : allStats) {
						stats.maxKMerCounts = countMap.get(stats.getTaxid());
					}
				}
			}
			this.uniqueCounter = null;
		} else {
			for (StatsPerTaxid stats : allStats) {
				stats.uniqueKmers = -1;
			}
		}

		return new Result(kmerStore.getK(), taxid2Stats, reads, kMers, countMap.get(null));
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
							for (int j = 0; j < entry.readSize; j++) {
								stats.maxContigDescriptor[j] = entry.readDescriptor[j];
							}
							stats.maxContigDescriptor[entry.readSize] = 0;
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
					synchronized (root) {
						stats = root.get(taxid, maxReadSize);
					}
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
							for (int j = 0; j < entry.readSize; j++) {
								stats.maxContigDescriptor[j] = entry.readDescriptor[j];
							}
							stats.maxContigDescriptor[entry.readSize] = 0;
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

		public MyReadEntry(int maxReadSizeBytes, boolean enablePrint) {
			super(maxReadSizeBytes);

			
			buffer = enablePrint ? new byte[maxReadSizeBytes * 4] : null; // It has to be rather long in some cases...
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

	protected static class TaxidStatsTrie extends DigitTrie<StatsPerTaxid> {
		@Override
		protected StatsPerTaxid createInGet(String digits, Object createContext) {
			return new StatsPerTaxid(digits, (int) createContext);
		}
	}

	public static class Result {
		private final int k;
		private final Map<String, StatsPerTaxid> taxid2Stats;
		private final long totalReads;
		private final long totalKMers;
		private final short[] totalMaxCounts;
		
		public Result(int k, Map<String, StatsPerTaxid> taxid2Stats, long totalReads, long totalKMers, short[] totalMaxCounts) {
			this.k = k;
			this.taxid2Stats = taxid2Stats;
			this.totalReads = totalReads;
			this.totalKMers = totalKMers;
			this.totalMaxCounts = totalMaxCounts;
		}
		
		public int getK() {
			return k;
		}

		public Map<String, StatsPerTaxid> getTaxid2Stats() {
			return Collections.unmodifiableMap(taxid2Stats);
		}

		public long getTotalKMers() {
			return totalKMers;
		}

		public long getTotalReads() {
			return totalReads;
		}
		
		public short[] getTotalMaxCounts() {
			return totalMaxCounts;
		}

		public boolean isWithMaxKMerCounts() {
			return totalMaxCounts != null;
		}
	}

	public static class StatsPerTaxid implements Serializable {
		private static final long serialVersionUID = 1L;
		
		protected String taxid;
		protected long reads;
		protected long uniqueKmers;
		protected long kmers;
		protected int maxContigLen;
		protected int contigs;
		protected short[] maxKMerCounts;
		protected byte[] maxContigDescriptor;

		public StatsPerTaxid(String taxid, int maxReadSizeBytes) {
			this.taxid = taxid;
			maxContigDescriptor = new byte[maxReadSizeBytes];
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
		
		public byte[] getMaxContigDescriptor() {
			return maxContigDescriptor;
		}

		public long getReads() {
			return reads;
		}

		public long getUniqueKMers() {
			return uniqueKmers;
		}

		public short[] getMaxKMerCounts() {
			return maxKMerCounts;
		}
	}
}
