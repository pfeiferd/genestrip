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

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.util.Arrays;
import java.util.Random;

import org.junit.Test;
import org.metagene.genestrip.DefaultExecutionContext;
import org.metagene.genestrip.ExecutionContext;
import org.metagene.genestrip.match.FastqKMerMatcher.MatcherReadEntry;
import org.metagene.genestrip.store.KMerStore;
import org.metagene.genestrip.store.KMerStore.ValueConverter;
import org.metagene.genestrip.store.RadixKMerStore;
import org.metagene.genestrip.tax.SmallTaxTree.SmallTaxIdNode;
import org.metagene.genestrip.util.CGAT;

public class HitCountingFastqKMerMatcherTest {
	private static final byte[] DECODE_TABLE = CGAT.newDecodeTable();
	private static final String[] TAXIDS = new String[] { "1", "2", "3" };
	private final Random random = new Random(42);

	// Counting in the database's entries.
	@Test
	public void testCountingInStore() {
		checkCounting(false);
	}

	// Counting in the matcher's own vector, as with parallelDbMatching=true.
	@Test
	public void testCountingInOwnVector() {
		checkCounting(true);
	}

	// Three 2-mers, one per tax id: CC (= GG), TT (= AA) and AG (= CT). Short random reads hit each of
	// them zero to several times, so the counts at the end must be min(hits, 3) per k-mer - and since
	// each tax id holds a single k-mer, a tax id's unique, once and twice counts follow directly.
	private void checkCounting(boolean own) {
		long ccKmer = CGAT.kMerToLong(new byte[] { 'C', 'C' }, 0, 2, null);
		long ttKmer = CGAT.kMerToLong(new byte[] { 'T', 'T' }, 0, 2, null);
		long agKmer = CGAT.kMerToLong(new byte[] { 'A', 'G' }, 0, 2, null);
		long[] kmers = new long[] { ccKmer, ttKmer, agKmer };
		int[] bucketSizes = new int[1 << RadixKMerStore.DEFAULT_RADIX_BITS];
		for (long kmer : kmers) {
			bucketSizes[RadixKMerStore.radixOf(kmer, RadixKMerStore.DEFAULT_RADIX_BITS)]++;
		}
		KMerStore<String> store = new RadixKMerStore<>(2, RadixKMerStore.DEFAULT_RADIX_BITS, bucketSizes, 0.0001,
				0.0001, Arrays.asList(TAXIDS), true);
		for (int i = 0; i < kmers.length; i++) {
			assertTrue(store.putLong(kmers[i], TAXIDS[i]));
		}
		store.optimize();

		ExecutionContext bundle = new DefaultExecutionContext(null, 0, 1000);
		int readLength = 6;
		CountingMatcher matcher = new CountingMatcher(store, readLength * 10, bundle, own);
		MatcherReadEntry entry = new MatcherReadEntry(2000, true, 4);
		entry.readSize = readLength;
		int[] hits = new int[TAXIDS.length];
		int[] seen = new int[4];

		for (int r = 0; r < 500; r++) {
			matcher.initStats();
			matcher.initUniqueCounter();
			Arrays.fill(hits, 0);
			// Several reads per run, so that hits on one k-mer accumulate across reads.
			for (int read = 0; read < 3; read++) {
				entry.bufferPos = 0;
				byte[] bases = entry.read;
				for (int j = 0; j < readLength; j++) {
					bases[j] = DECODE_TABLE[random.nextInt(4)];
					if (j > 0) {
						byte a = bases[j - 1];
						byte b = bases[j];
						if ((a == 'C' && b == 'C') || (a == 'G' && b == 'G')) {
							hits[0]++;
						} else if ((a == 'T' && b == 'T') || (a == 'A' && b == 'A')) {
							hits[1]++;
						} else if ((a == 'A' && b == 'G') || (a == 'C' && b == 'T')) {
							hits[2]++;
						}
					}
				}
				matcher.matchRead(entry, 0);
			}
			matcher.computeUniqueKmerCounts();
			for (int j = 0; j < TAXIDS.length; j++) {
				CountsPerTaxid stats = matcher.getStats(TAXIDS[j]);
				if (hits[j] == 0) {
					assertNull(stats);
				} else {
					seen[Math.min(hits[j], 3)]++;
					assertEquals(hits[j], stats.getKMers());
					assertEquals(1, stats.getUniqueKMers());
					assertEquals(hits[j] == 1 ? 1L : 0L, (long) stats.getKMersHitOnce());
					assertEquals(hits[j] == 2 ? 1L : 0L, (long) stats.getKMersHitTwice());
				}
			}
		}
		// The reads must have produced all three non-zero counter states.
		for (int c = 1; c <= 3; c++) {
			assertTrue("counter state " + c + " never reached", seen[c] > 0);
		}
		assertEquals(!own, store.isMarkVisited());
		matcher.dump();
	}

	private static class CountingMatcher extends HitCountingFastqKMerMatcher {
		private final KMerStore<String> orgKmerStore;
		private final boolean own;

		public CountingMatcher(final KMerStore<String> kmerStore, int initialReadSize, ExecutionContext bundle,
				boolean own) {
			super(kmerStore.convertValues(new ValueConverter<String, SmallTaxIdNode>() {
				@Override
				public SmallTaxIdNode convertValue(String value) {
					SmallTaxIdNode node = new SmallTaxIdNode(value, null, null);
					node.setStoreIndex(kmerStore.getIndexForValue(value));
					return node;
				}
			}), initialReadSize, 1, bundle, false, null, 4, 0, 0, true, 1, null);
			orgKmerStore = kmerStore;
			this.own = own;
		}

		@Override
		protected boolean isOwnUniqueKMerBits() {
			return own;
		}

		public CountsPerTaxid getStats(String taxid) {
			return statsIndex[orgKmerStore.getIndexForValue(taxid)];
		}
	}
}
