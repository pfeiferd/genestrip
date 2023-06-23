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
package org.metagene.genestrip.trie;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.codec.digest.MurmurHash3;
import org.metagene.genestrip.fastq.AbstractFastqReader;
import org.metagene.genestrip.util.ByteArrayUtil;
import org.metagene.genestrip.util.StreamProvider;
import org.metagene.genestrip.util.StreamProvider.ByteCountingInputStreamAccess;

import it.unimi.dsi.fastutil.ints.Int2ObjectLinkedOpenHashMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectMap;
import it.unimi.dsi.fastutil.objects.ObjectCollection;

public class FastqTrieClassifier extends AbstractFastqReader {
	private final KMerTrie<String> trie;
	private final int k;
	private final KmerDuplicationCount duplicationCount;

	private TaxidStatsTrie root;

	private ByteCountingInputStreamAccess byteCountAccess;
	private long fastqFileSize;
	private long startTime;
	private long indexedC;

	private OutputStream indexed;

	// A PrintStream is implicitly synchronized. So we don't need to worry about
	// multi threading when using it.
	private PrintStream out;

	public FastqTrieClassifier(KMerTrie<String> trie, int maxReadSize, int maxQueueSize, int consumerNumber,
			boolean withDupCount) {
		super(maxReadSize, maxQueueSize, consumerNumber);
		this.trie = trie;
		this.k = trie.getLen();
		duplicationCount = withDupCount ? new KmerDuplicationCount(k, 42) : null;
	}

	@Override
	protected ReadEntry createReadEntry(int maxReadSizeBytes) {
		return new MyReadEntry(maxReadSizeBytes);
	}

	public List<StatsPerTaxid> runClassifier(File fastq, File filteredFile, File krakenOutStyleFile)
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
		indexedC = 0;

		root = new TaxidStatsTrie();
		if (duplicationCount != null) {
			duplicationCount.clear();
		}

		readFastq(byteCountAccess.getInputStream());

		if (indexed != null) {
			indexed.close();
			out.close();
		}
		byteCountAccess.getInputStream().close();

		List<StatsPerTaxid> counts = new ArrayList<FastqTrieClassifier.StatsPerTaxid>();
		root.collectValues(counts);

		if (duplicationCount != null) {
			for (StatsPerTaxid stats : counts) {
				stats.uniqueKmers = duplicationCount.getUniqeKmers(stats.taxid);
				if (stats.kmers != duplicationCount.getKmerCount(stats.taxid)) {
					if (logger.isInfoEnabled()) {
						logger.info("Warning: inconsistent kmer counts for taxid: " + stats.taxid);
					}
				}
			}
		} else {
			for (StatsPerTaxid stats : counts) {
				stats.uniqueKmers = -1;
			}
		}
		return counts;
	}

	@Override
	protected void nextEntry(ReadEntry entry) throws IOException {
		MyReadEntry myEntry = (MyReadEntry) entry;
		myEntry.bufferPos = 0;

		boolean found = classifyRead(myEntry, false);
		if (!found) {
			found = classifyRead(myEntry, true);
		}
		if (found) {
			if (indexed != null) {
				rewriteInput(entry, indexed);
			}
			if (out != null) {
				myEntry.printChar('\n');
				synchronized (out) {
					out.write(myEntry.buffer, 0, myEntry.bufferPos);
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
		if (logger.isInfoEnabled()) {
			if (reads % 100000 == 0) {
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

//	protected boolean classifyRead(ReadEntry entry, boolean reverse) {
//		String res = null;
//		int max = entry.readSize - k;
//		boolean found = false;
//		for (int i = 0; i <= max; i++) {
//			res = trie.get(entry.read, i, reverse);
//			if (res != null) {
//				found = true;
//				root.inc(res);
//			}
//		}
//		return found;
//	}

	protected boolean classifyRead(MyReadEntry entry, boolean reverse) {
		boolean found = false;
		int prints = 0;

		String taxid = null;
		int max = entry.readSize - k;
		String lastTaxid = null;
		int contigLen = 0;
		StatsPerTaxid stats = null;

		for (int i = 0; i <= max; i++) {
			taxid = trie.get(entry.read, i, reverse);
			if (contigLen > 0 && taxid != lastTaxid) {
				if (out != null) {
					printKrakenStyleOut(entry, lastTaxid, contigLen, prints++, reverse);
				}
				contigLen = 0;
			}
			if (taxid != null) {
				stats = root.get(taxid);
				if (stats == null) {
					stats = root.create(taxid);
				}
				synchronized (stats) {
					if (!found) {
						stats.reads++;
					}
					found = true;
					stats.kmers++;
					if (taxid != lastTaxid) {
						stats.contigs++;
						stats.sumContigsLen += contigLen;
						if (contigLen > stats.maxContigLen) {
							stats.maxContigLen = contigLen;
						}
					}
				}
				if (duplicationCount != null) {
					if (!duplicationCount.put(taxid, entry.read, i, reverse)) {
						if (logger.isInfoEnabled()) {
							synchronized (logger) {
								logger.info("Warning: Updating duplication count failed for kmer under taxid " + taxid);
							}
						}
					}
				}
			}
			contigLen++;
			lastTaxid = taxid;
		}
		if (found && contigLen > 0) {
			if (taxid != null) {
				stats.contigs++;
				stats.sumContigsLen += contigLen;
				if (contigLen > stats.maxContigLen) {
					stats.maxContigLen = contigLen;
				}
			}
			if (out != null) {
				printKrakenStyleOut(entry, lastTaxid, contigLen, prints, reverse);
			}
		}

		return found;
	}

	protected void printKrakenStyleOut(MyReadEntry entry, String taxid, int contigLen, int state, boolean reverse) {
		if (state == 0) {
			entry.printBytes(entry.readDescriptor);
			entry.printChar('\t');
			entry.printChar('0');
			entry.printChar('\t');
			entry.printInt(entry.readSize);
			entry.printChar('\t');
			// entry.printChar(reverse ? 'R' : 'S');
		} else {
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

	protected static class MyReadEntry extends ReadEntry {
		public byte[] buffer;
		public int bufferPos;

		public MyReadEntry(int maxReadSizeBytes) {
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
	}

	public static class KmerDuplicationCount {
		private final Int2ObjectMap<Entry> map;
		private final int k;
		private final int seed;

		public KmerDuplicationCount(int k, int seed) {
			this.k = k;
			this.seed = seed;
			map = new Int2ObjectLinkedOpenHashMap<Entry>();
		}

		public void clear() {
			map.clear();
		}

		protected boolean put(String taxid, byte[] read, int start, boolean reverse) {
			int hash = MurmurHash3.hash32x86(read, start, k, seed); // FIX me this does not respect reverse reads!

			synchronized (map) {
				Entry e = map.get(hash);
				if (e == null) {
					e = new Entry();
					e.taxid = taxid;
					map.put(hash, e);
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

		public int getKmerCount(String taxid) {
			int res = 0;
			ObjectCollection<Entry> entries = map.values();
			for (Entry e : entries) {
				if (e.taxid == taxid) {
					res += e.count;
				}
			}
			return res;
		}

		private static class Entry {
			public int count;
			public String taxid;
		}
	}

	public static class StatsPerTaxid {
		protected String taxid;
		protected int reads;
		protected int uniqueKmers;
		protected int kmers;
		protected int sumContigsLen;
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

		public int getKmers() {
			return kmers;
		}

		public int getMaxContigLen() {
			return maxContigLen;
		}

		public int getReads() {
			return reads;
		}

		public int getSumContigsLen() {
			return sumContigsLen;
		}

		public int getUniqueKmers() {
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
}
