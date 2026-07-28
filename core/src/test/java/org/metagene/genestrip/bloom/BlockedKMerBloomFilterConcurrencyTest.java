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
package org.metagene.genestrip.bloom;

import static org.junit.Assert.assertTrue;

import java.util.Random;
import java.util.concurrent.CyclicBarrier;

import org.junit.Test;

/**
 * Verifies that {@link BlockedKMerBloomFilter#putLong(long)} is safe for concurrent inserts on both
 * backings: filling a filter from several threads must not lose any bits, i.e. every inserted key is
 * still found afterwards. Without locking the {@code word |= mask} read-modify-write can drop a
 * concurrent update to the same word, which shows up as a false negative here.
 */
public class BlockedKMerBloomFilterConcurrencyTest {

	private static final int THREADS = 8;
	private static final int KMERS_PER_THREAD = 20000;

	@Test
	public void testConcurrentPutSmallBacking() throws Exception {
		assertNoFalseNegativesAfterConcurrentFill(newFilter(false));
	}

	@Test
	public void testConcurrentPutLargeBacking() throws Exception {
		assertNoFalseNegativesAfterConcurrentFill(newFilter(true));
	}

	/**
	 * Creates a filter sized for the whole key set, with a bucket width small enough that keys spread
	 * over many buckets (and so over many locks) when the bucketed backing is used.
	 */
	private static BlockedKMerBloomFilter newFilter(boolean large) {
		long size = (long) THREADS * KMERS_PER_THREAD;
		// The filter picks its backing from the size relative to MAX_SMALL_CAPACITY at construction,
		// which this test size never reaches, so ask for the bucketed backing explicitly.
		int bucketShift = BlockedKMerBloomFilter.MIN_BUCKET_SHIFT + 4;
		return large ? BlockedKMerBloomFilter.newLargeBacked(size, 10, 42L, bucketShift)
				: new BlockedKMerBloomFilter(size, 10, 42L, bucketShift);
	}

	/**
	 * Fills the filter from {@link #THREADS} threads that start together, each inserting its own keys,
	 * and asserts that every key is found once all of them have finished.
	 */
	private static void assertNoFalseNegativesAfterConcurrentFill(BlockedKMerBloomFilter filter) throws Exception {
		long[][] kmers = new long[THREADS][];
		for (int t = 0; t < THREADS; t++) {
			kmers[t] = randomKMers(KMERS_PER_THREAD, t);
		}
		// The barrier makes the threads insert at the same time, so their words actually contend.
		CyclicBarrier start = new CyclicBarrier(THREADS);
		Thread[] threads = new Thread[THREADS];
		for (int t = 0; t < THREADS; t++) {
			long[] mine = kmers[t];
			threads[t] = new Thread(() -> {
				try {
					start.await();
				} catch (Exception e) {
					throw new IllegalStateException(e);
				}
				for (long kmer : mine) {
					filter.putLong(kmer);
				}
			});
			threads[t].start();
		}
		for (Thread thread : threads) {
			thread.join();
		}
		// Joining every writer establishes happens-before, so these reads see all inserts.
		for (long[] mine : kmers) {
			for (long kmer : mine) {
				assertTrue("inserted key must be found", filter.containsLong(kmer));
			}
		}
	}

	private static long[] randomKMers(int n, long seed) {
		Random random = new Random(seed);
		long[] kmers = new long[n];
		for (int i = 0; i < n; i++) {
			kmers[i] = random.nextLong();
		}
		return kmers;
	}
}
