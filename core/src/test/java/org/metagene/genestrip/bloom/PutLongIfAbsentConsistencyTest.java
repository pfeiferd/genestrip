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

import static org.junit.Assert.assertArrayEquals;
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
import org.metagene.genestrip.util.LargeBitVector;

/**
 * Verifies that the lock-free combined insert {@link AbstractKMerBloomFilter#putLongIfAbsent(long)}
 * yields a filter that is bit-for-bit identical to the classic {@code if (!containsLong(x))
 * putLong(x)} sequence, both single-threaded (small and large backing) and under concurrent
 * insertion.
 */
public class PutLongIfAbsentConsistencyTest {
	private static final double FPP = 0.001;

	private static XORKMerBloomFilter newFilter(long expected, boolean large) {
		XORKMerBloomFilter filter = new XORKMerBloomFilter(FPP);
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

	private static void assertBitVectorsEqual(LargeBitVector expected, LargeBitVector actual) {
		assertEquals("bit size", expected.getBitSize(), actual.getBitSize());
		assertEquals("large flag", expected.isLarge(), actual.isLarge());
		long bitSize = expected.getBitSize();
		for (long b = 0; b < bitSize; b++) {
			if (expected.get(b) != actual.get(b)) {
				throw new AssertionError("bit " + b + " differs");
			}
		}
		// Extra exactness check for the small backing where the whole word array is directly available.
		if (expected.bits != null) {
			assertArrayEquals(expected.bits, actual.bits);
		}
	}

	private void checkMatchesClassicPath(boolean large) {
		int n = 20_000;
		long[] kmers = randomKMers(n, 12345);

		XORKMerBloomFilter reference = newFilter(n, large);
		for (long kmer : kmers) {
			// Feed each k-mer twice to exercise the "already present" branch.
			if (!reference.containsLong(kmer)) {
				reference.putLong(kmer);
			}
			if (!reference.containsLong(kmer)) {
				reference.putLong(kmer);
			}
		}

		XORKMerBloomFilter candidate = newFilter(n, large);
		for (long kmer : kmers) {
			boolean firstAdded = candidate.putLongIfAbsent(kmer);
			boolean secondAdded = candidate.putLongIfAbsent(kmer);
			// The first insert of a genuinely new k-mer is reported as added; the immediate repeat is
			// never reported as added.
			assertFalse("repeat insert must not be reported as new", secondAdded);
			if (!firstAdded) {
				// If not reported new it must have already been present (a hash collision with an
				// earlier k-mer); membership must then hold.
				assertTrue(candidate.containsLong(kmer));
			}
		}

		assertBitVectorsEqual(reference.bitVector, candidate.bitVector);
		assertEquals("entries", reference.getEntries(), candidate.getEntries());
		for (long kmer : kmers) {
			assertTrue("no false negative", candidate.containsLong(kmer));
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

	@Test
	public void testSerializationFoldsConcurrentEntries() throws IOException, ClassNotFoundException {
		int n = 5_000;
		long[] kmers = randomKMers(n, 4242);
		XORKMerBloomFilter filter = newFilter(n, false);
		for (long kmer : kmers) {
			filter.putLongIfAbsent(kmer);
		}
		long entriesBefore = filter.getEntries();
		assertTrue("some entries expected", entriesBefore > 0);

		ByteArrayOutputStream bytes = new ByteArrayOutputStream();
		try (ObjectOutputStream out = new ObjectOutputStream(bytes)) {
			out.writeObject(filter);
		}
		XORKMerBloomFilter loaded;
		try (ObjectInputStream in = new ObjectInputStream(new ByteArrayInputStream(bytes.toByteArray()))) {
			loaded = (XORKMerBloomFilter) in.readObject();
		}

		// The concurrent counter is transient, so its value must have been folded into 'entries' on
		// write for the reloaded filter to report the same count.
		assertEquals("entries survive serialization", entriesBefore, loaded.getEntries());
		for (long kmer : kmers) {
			assertTrue("no false negative after reload", loaded.containsLong(kmer));
		}
	}

	@Test
	public void testConcurrentInsertMatchesSingleThreaded() throws InterruptedException {
		int threadCount = 8;
		int perThread = 25_000;
		int total = threadCount * perThread;
		// Distinct k-mers partitioned across threads (no key inserted by two threads), so the "newly
		// added" count is exact even under the benign concurrent-overcount race.
		long[] kmers = randomKMers(total, 99);

		XORKMerBloomFilter reference = newFilter(total, false);
		for (long kmer : kmers) {
			reference.putLongIfAbsent(kmer);
		}

		XORKMerBloomFilter concurrent = newFilter(total, false);
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

		// Order-independent OR of bits: concurrent build must be bit-identical to the serial build,
		// which proves the two filters answer every membership query the same way.
		assertBitVectorsEqual(reference.bitVector, concurrent.bitVector);
		for (long kmer : kmers) {
			assertTrue("no false negative after concurrent insert", concurrent.containsLong(kmer));
		}
		// The distinct k-mers are partitioned across threads, so essentially all are counted as new;
		// the only slack is the rare k-mer whose bits are fully covered by collisions, whose "new"
		// verdict is insertion-order sensitive. Bound it loosely rather than requiring exact equality.
		long entries = concurrent.getEntries();
		assertTrue("entries too low: " + entries, entries >= total - total / 100);
		assertTrue("entries above insert count: " + entries, entries <= total);
	}
}
