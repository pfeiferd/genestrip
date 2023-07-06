package org.metagene.genestrip.trie;

import static org.junit.Assert.assertEquals;

import java.util.Arrays;
import java.util.Map;
import java.util.Random;

import org.junit.Test;
import org.metagene.genestrip.store.FastqKMerMatcher;
import org.metagene.genestrip.store.FastqKMerMatcher.MyReadEntry;
import org.metagene.genestrip.store.FastqKMerMatcher.StatsPerTaxid;
import org.metagene.genestrip.store.KMerStore;
import org.metagene.genestrip.util.ByteArrayUtil;
import org.metagene.genestrip.util.CGAT;
import org.metagene.genestrip.util.CGATRingBuffer;

public class FastqMatcherTest {
	private Random random = new Random(42);
	private static final String[] taxids = new String[] { "1", "2", "3" };

	@Test
	public void testMatchRead() {
		int readLength = 150;
		int entries = 200;

		MyFastqMatcher matcher = new MyFastqMatcher(new SimpleTestKmerStore(), readLength, 1, 0, true);

		MyReadEntry entry = new MyReadEntry(2000, 1);
		entry.readSize = readLength;

		long[] counters = new long[taxids.length];
		int[] contigs = new int[taxids.length];
		int[] maxContigLen = new int[taxids.length];
		int previousPos = 0;
		int contigLen;

		for (int i = 1; i <= entries; i++) {
			// We check correctness of stats for each read separately.
			Arrays.fill(counters, 0);
			Arrays.fill(contigs, 0);
			Arrays.fill(maxContigLen, 0);
			matcher.initRoot();
			matcher.getDuplicationCount().clear();
			entry.bufferPos = 0;
			contigLen = 0;

			for (int j = 0; j < readLength; j++) {
				int pos = random.nextInt(4);
				entry.read[j] = CGAT.DECODE_TABLE[pos];
				System.out.print((char) entry.read[j]);
				if (pos != 3) {
					counters[pos]++;
				}
				if (j > 0) {
					if (pos != previousPos) {
						if (previousPos != 3) {
							contigs[previousPos]++;
							if (contigLen > maxContigLen[previousPos]) {
								maxContigLen[previousPos] = contigLen;
							}
						}
						contigLen = 0;
					}
				}
				contigLen++;
				previousPos = pos;
			}
			if (previousPos != 3) {
				contigs[previousPos]++;
				if (contigLen > maxContigLen[previousPos]) {
					maxContigLen[previousPos] = contigLen;
				}
			}
			ByteArrayUtil.print(entry.read, System.out);
			matcher.classifyRead(entry, false);

			ByteArrayUtil.print(entry.read, System.out);
			System.out.println();
			System.out.write(entry.buffer, 0, entry.bufferPos);
			System.out.println();

			for (int j = 0; j < taxids.length; j++) {
				StatsPerTaxid stats = matcher.getStats(taxids[j]);
				assertEquals(stats.getKMers(), counters[j]);
				assertEquals(counters[j], matcher.getDuplicationCount().getKmerCount((long) j));
				assertEquals(contigs[j], stats.getContigs());
				assertEquals(maxContigLen[j], stats.getMaxContigLen());
			}
		}
	}

	protected static class SimpleTestKmerStore implements KMerStore<String> {
		private static final long serialVersionUID = 1L;

		@Override
		public int getMaxValues() {
			return -1;
		}

		@Override
		public int getK() {
			return 1;
		}

		@Override
		public long getEntries() {
			return 4;
		}

		@Override
		public void initSize(long size) {
		}
				
		@Override
		public Map<String, Long> getNKmersPerTaxid() {
			throw new UnsupportedOperationException("Not implemented for this test.");
		}
		

		@Override
		public long getSize() {
			return Long.MAX_VALUE;
		}

		@Override
		public boolean put(CGATRingBuffer buffer, String value, boolean reverse) {
			throw new UnsupportedOperationException();
		}

		@Override
		public boolean put(byte[] nseq, int start, String value, boolean reverse) {
			throw new UnsupportedOperationException();
		}

		@Override
		public String get(CGATRingBuffer buffer, boolean reverse) {
			if (buffer.get(0) == 'T') {
				return null;
			}
			return taxids[CGAT.CGAT_JUMP_TABLE[buffer.get(0)]];
		}

		@Override
		public String get(byte[] nseq, int start, boolean reverse) {
			if (nseq[start] == 'T') {
				return null;
			}
			return taxids[CGAT.CGAT_JUMP_TABLE[nseq[start]]];
		}

		@Override
		public void visit(KMerStoreVisitor<String> visitor) {
			throw new UnsupportedOperationException();
		}

		@Override
		public void optimize() {
		}

		@Override
		public boolean isOptimized() {
			return true;
		}
	}

	protected static class MyFastqMatcher extends FastqKMerMatcher {
		public MyFastqMatcher(KMerStore<String> kmerStore, int maxReadSize, int maxQueueSize, int consumerNumber,
				boolean withDupCount) {
			super(kmerStore, maxReadSize, maxQueueSize, consumerNumber, withDupCount);
			out = System.out;
		}

		@Override
		public void initRoot() {
			super.initRoot();
		}

		public StatsPerTaxid getStats(String taxid) {
			return root.get(taxid);
		}

		@Override
		public boolean classifyRead(MyReadEntry entry, boolean reverse) {
			return super.classifyRead(entry, reverse);
		}

		public KmerDuplicationCount getDuplicationCount() {
			return duplicationCount;
		}
	}
}
