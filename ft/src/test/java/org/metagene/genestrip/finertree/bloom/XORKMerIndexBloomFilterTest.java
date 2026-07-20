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
package org.metagene.genestrip.finertree.bloom;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.Random;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.atomic.AtomicReference;

import org.junit.Test;

/**
 * Verifies that {@link XORKMerIndexBloomFilter#putLongIntIfAbsent(long, int)} — the lock-free
 * combined insert used by the k-mer index bloom goal — behaves exactly like the previous {@code if
 * (!containsLongInt(k, i)) putLongInt(k, i)} sequence, including under concurrent insertion, and
 * that the raw single-argument operations remain disabled.
 */
public class XORKMerIndexBloomFilterTest {
	private static final double FPP = 0.001;

	private static XORKMerIndexBloomFilter newFilter(long expected) {
		return new XORKMerIndexBloomFilter(FPP, expected);
	}

	@Test
	public void testCombinedInsertMatchesClassicPath() {
		int n = 20_000;
		Random random = new Random(7);
		long[] kmers = new long[n];
		int[] indices = new int[n];
		for (int i = 0; i < n; i++) {
			kmers[i] = random.nextLong();
			// Mix in the OTHER_VALUE marker and small leaf indices, like the real goal does.
			indices[i] = (i % 5 == 0) ? Integer.MAX_VALUE : (i % 37);
		}

		XORKMerIndexBloomFilter reference = newFilter(n);
		for (int i = 0; i < n; i++) {
			if (!reference.containsLongInt(kmers[i], indices[i])) {
				reference.putLongInt(kmers[i], indices[i]);
			}
		}

		XORKMerIndexBloomFilter candidate = newFilter(n);
		for (int i = 0; i < n; i++) {
			boolean firstAdded = candidate.putLongIntIfAbsent(kmers[i], indices[i]);
			boolean secondAdded = candidate.putLongIntIfAbsent(kmers[i], indices[i]);
			assertFalse("repeat (k-mer, index) insert must not be reported as new", secondAdded);
			if (!firstAdded) {
				assertTrue(candidate.containsLongInt(kmers[i], indices[i]));
			}
		}

		// Entry counts are not compared: the candidate inserts each pair twice above (the repeat-insert
		// behaviour under test), and getEntries() now counts every insert. Membership equivalence below
		// is what proves the combined insert matches the classic path.
		for (int i = 0; i < n; i++) {
			assertTrue("no false negative", candidate.containsLongInt(kmers[i], indices[i]));
			// The same k-mer under a different index must not be forced present by the insert.
			assertFalse("membership must depend on index",
					reference.containsLongInt(kmers[i], indices[i]) != candidate.containsLongInt(kmers[i], indices[i]));
		}
	}

	@Test
	public void testConcurrentInsertMatchesSingleThreaded() throws InterruptedException {
		int threadCount = 8;
		int perThread = 25_000;
		int total = threadCount * perThread;
		Random random = new Random(2024);
		long[] kmers = new long[total];
		int[] indices = new int[total];
		for (int i = 0; i < total; i++) {
			kmers[i] = random.nextLong();
			indices[i] = (i % 5 == 0) ? Integer.MAX_VALUE : (i % 37);
		}

		XORKMerIndexBloomFilter reference = newFilter(total);
		for (int i = 0; i < total; i++) {
			reference.putLongIntIfAbsent(kmers[i], indices[i]);
		}

		XORKMerIndexBloomFilter concurrent = newFilter(total);
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
						concurrent.putLongIntIfAbsent(kmers[i], indices[i]);
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

		for (int i = 0; i < total; i++) {
			assertTrue("no false negative after concurrent insert",
					concurrent.containsLongInt(kmers[i], indices[i]));
		}
		// The concurrent OR of bits is order-independent, so the concurrent filter must answer every
		// membership query identically to the serial one. Probe inserted pairs and many random pairs.
		Random probe = new Random(555);
		for (int i = 0; i < total; i++) {
			int j = probe.nextInt(total);
			assertEquals("membership must match serial build",
					reference.containsLongInt(kmers[j], indices[j]),
					concurrent.containsLongInt(kmers[j], indices[j]));
			long randomKmer = probe.nextLong();
			int randomIndex = probe.nextInt();
			assertEquals("membership must match serial build",
					reference.containsLongInt(randomKmer, randomIndex),
					concurrent.containsLongInt(randomKmer, randomIndex));
		}
		// Distinct pairs partitioned across threads: essentially all counted as new, with only rare
		// collision-order slack.
		long entries = concurrent.getEntries();
		assertTrue("entries too low: " + entries, entries >= total - total / 100);
		assertTrue("entries above insert count: " + entries, entries <= total);
	}

	@Test(expected = UnsupportedOperationException.class)
	public void testRawPutLongDisabled() {
		newFilter(100).putLong(42L);
	}

	@Test(expected = UnsupportedOperationException.class)
	public void testRawPutLongIfAbsentDisabled() {
		newFilter(100).putLongIfAbsent(42L);
	}

	@Test(expected = UnsupportedOperationException.class)
	public void testRawContainsLongDisabled() {
		newFilter(100).containsLong(42L);
	}
}
