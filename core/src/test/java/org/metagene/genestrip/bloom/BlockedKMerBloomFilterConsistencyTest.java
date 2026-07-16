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

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.util.Random;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.atomic.AtomicReference;

import org.junit.Test;

/**
 * Verifies that {@link BlockedKMerBloomFilter#putLongIfAbsent(long)} — the lock-free combined insert
 * used by the temporary-index size estimator — behaves exactly like the previous {@code if
 * (!containsLong(x)) putLong(x)} sequence, single-threaded and under concurrent insertion, and that
 * the entry count survives serialization.
 */
public class BlockedKMerBloomFilterConsistencyTest {

	private static BlockedKMerBloomFilter newFilter(long expected) {
		return newFilter(expected, false);
	}

	private static BlockedKMerBloomFilter newFilter(long expected, boolean large) {
		BlockedKMerBloomFilter filter = new BlockedKMerBloomFilter();
		filter.ensureExpectedSize(expected, large);
		return filter;
	}

	private static long[] randomKMers(int n, long seed) {
		Random random = new Random(seed);
		long[] kmers = new long[n];
		for (int i = 0; i < n; i++) {
			kmers[i] = random.nextLong();
		}
		return kmers;
	}

	/** Asserts the two filters answer every membership query identically over a large probe set. */
	private static void assertEquivalent(BlockedKMerBloomFilter expected, BlockedKMerBloomFilter actual,
			long[] insertedKeys, long probeSeed) {
		for (long kmer : insertedKeys) {
			assertEquals("membership on inserted key", expected.containsLong(kmer), actual.containsLong(kmer));
		}
		Random probe = new Random(probeSeed);
		for (int i = 0; i < insertedKeys.length; i++) {
			long randomKmer = probe.nextLong();
			assertEquals("membership on random key", expected.containsLong(randomKmer), actual.containsLong(randomKmer));
		}
	}

	@Test
	public void testMatchesClassicPathSmallBacking() {
		checkMatchesClassicPath(false);
	}

	@Test
	public void testMatchesClassicPathLargeBacking() {
		checkMatchesClassicPath(true);
	}

	private void checkMatchesClassicPath(boolean large) {
		int n = 20_000;
		long[] kmers = randomKMers(n, 314);

		BlockedKMerBloomFilter reference = newFilter(n, large);
		for (long kmer : kmers) {
			// Feed each k-mer twice to exercise the "already present" branch.
			if (!reference.containsLong(kmer)) {
				reference.putLong(kmer);
			}
			if (!reference.containsLong(kmer)) {
				reference.putLong(kmer);
			}
		}

		BlockedKMerBloomFilter candidate = newFilter(n, large);
		for (long kmer : kmers) {
			boolean firstAdded = candidate.putLongIfAbsent(kmer);
			boolean secondAdded = candidate.putLongIfAbsent(kmer);
			assertFalse("repeat insert must not be reported as new", secondAdded);
			if (!firstAdded) {
				assertTrue(candidate.containsLong(kmer));
			}
		}

		// The two backing words a key touches are always distinct, so single-threaded "newly added"
		// detection is exact and both the entry count and every membership answer must match.
		assertEquals("entries", reference.getEntries(), candidate.getEntries());
		assertEquivalent(reference, candidate, kmers, 999);
	}

	@Test
	public void testConcurrentInsertMatchesSingleThreaded() throws InterruptedException {
		int threadCount = 8;
		int perThread = 25_000;
		int total = threadCount * perThread;
		long[] kmers = randomKMers(total, 27);

		BlockedKMerBloomFilter reference = newFilter(total);
		for (long kmer : kmers) {
			reference.putLongIfAbsent(kmer);
		}

		BlockedKMerBloomFilter concurrent = newFilter(total);
		CountDownLatch start = new CountDownLatch(1);
		Thread[] threads = new Thread[threadCount];
		AtomicReference<Throwable> failure = new AtomicReference<>();
		for (int t = 0; t < threadCount; t++) {
			final int from = t * perThread;
			final int to = from + perThread;
			threads[t] = new Thread(() -> {
				try {
					start.await();
					for (int i = from; i < to; i++) {
						concurrent.putLongIfAbsent(kmers[i]);
					}
				} catch (Throwable th) {
					failure.compareAndSet(null, th);
				}
			});
			threads[t].start();
		}
		start.countDown();
		for (Thread thread : threads) {
			thread.join();
		}
		assertEquals(null, failure.get());

		// The concurrent OR of bits is order-independent, so both filters must answer identically.
		assertEquivalent(reference, concurrent, kmers, 555);
		// Distinct keys partitioned across threads: essentially all counted as new, with only rare
		// collision-order slack.
		long entries = concurrent.getEntries();
		assertTrue("entries too low: " + entries, entries >= total - total / 100);
		assertTrue("entries above insert count: " + entries, entries <= total);
	}

	@Test
	public void testSerializationFoldsConcurrentEntries() throws IOException, ClassNotFoundException {
		int n = 5_000;
		long[] kmers = randomKMers(n, 4242);
		BlockedKMerBloomFilter filter = newFilter(n);
		for (long kmer : kmers) {
			filter.putLongIfAbsent(kmer);
		}
		long entriesBefore = filter.getEntries();
		assertTrue("some entries expected", entriesBefore > 0);

		ByteArrayOutputStream bytes = new ByteArrayOutputStream();
		try (ObjectOutputStream out = new ObjectOutputStream(bytes)) {
			out.writeObject(filter);
		}
		BlockedKMerBloomFilter loaded;
		try (ObjectInputStream in = new ObjectInputStream(new ByteArrayInputStream(bytes.toByteArray()))) {
			loaded = (BlockedKMerBloomFilter) in.readObject();
		}

		assertEquals("entries survive serialization", entriesBefore, loaded.getEntries());
		for (long kmer : kmers) {
			assertTrue("no false negative after reload", loaded.containsLong(kmer));
		}
	}
}
