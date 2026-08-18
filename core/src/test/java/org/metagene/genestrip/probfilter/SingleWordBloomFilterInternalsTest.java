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
package org.metagene.genestrip.probfilter;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.concurrent.CountDownLatch;

import org.junit.Test;

/**
 * Tests what is specific to {@link SingleWordBloomFilter} beyond the shared
 * {@link BloomFilterTest} suite: that a key really only ever touches one word, that its bits
 * survive whatever else is inserted, and the sizing, serialization and concurrency behaviour.
 */
public class SingleWordBloomFilterInternalsTest {
	private static final int SIZE = 200000;

	private long[] randomKeys(int n, long seed) {
		Random random = new Random(seed);
		long[] keys = new long[n];
		for (int i = 0; i < n; i++) {
			keys[i] = random.nextLong();
		}
		return keys;
	}

	@Test
	public void testNoFalseNegativesAcrossSizings() {
		// The contract every ProbFilter owes, checked over the sizings and both backings. Bits are
		// only ever ORed in, so no insert can ever remove another key's bits.
		for (int bitsPerKey : new int[] { 8, 10, 16, 24 }) {
			for (boolean large : new boolean[] { false, true }) {
				SingleWordBloomFilter filter = large ? SingleWordBloomFilter.newLargeBacked(SIZE, bitsPerKey)
						: new SingleWordBloomFilter(SIZE, bitsPerKey);
				assertEquals("large backing", large, filter.isLargeBacked());
				long[] keys = randomKeys(SIZE, 42);
				for (long key : keys) {
					filter.putLong(key);
				}
				for (long key : keys) {
					assertTrue("bitsPerKey=" + bitsPerKey + " large=" + large, filter.containsLong(key));
				}
			}
		}
	}

	@Test
	public void testKeyTouchesExactlyOneWord() {
		// The point of the design: inserting one key must leave every word but its own untouched.
		for (int hashBits : new int[] { 1, 5, SingleWordBloomFilter.MAX_HASH_BITS }) {
			SingleWordBloomFilter filter = new SingleWordBloomFilter(700, 64, hashBits);
			long words = filter.getBitSize() / 64;
			assertTrue("expected several words, got " + words, words > 1);
			for (long key : randomKeys(200, 17)) {
				filter.clear();
				filter.putLong(key);
				int touched = 0;
				for (long index = 0; index < words; index++) {
					if (filter.getWord(index) != 0) {
						touched++;
					}
				}
				assertEquals("hashBits=" + hashBits, 1, touched);
			}
		}
	}

	@Test
	public void testKeySetsAtMostHashBitsBits() {
		// A key sets hashBits bits, or fewer when two of its bit positions coincide - never more.
		for (int hashBits : new int[] { 1, 3, 6, SingleWordBloomFilter.MAX_HASH_BITS }) {
			SingleWordBloomFilter filter = new SingleWordBloomFilter(700, 64, hashBits);
			long words = filter.getBitSize() / 64;
			for (long key : randomKeys(200, 23)) {
				filter.clear();
				filter.putLong(key);
				int bits = 0;
				for (long index = 0; index < words; index++) {
					bits += Long.bitCount(filter.getWord(index));
				}
				assertTrue("hashBits=" + hashBits + " set " + bits + " bits", bits >= 1 && bits <= hashBits);
			}
		}
	}

	@Test
	public void testPutLongReportsNewlyAdded() {
		SingleWordBloomFilter filter = new SingleWordBloomFilter(SIZE, 24);
		long[] keys = randomKeys(1000, 4711);
		for (long key : keys) {
			assertTrue(filter.putLong(key));
		}
		for (long key : keys) {
			assertFalse(filter.putLong(key));
		}
	}

	@Test
	public void testOptimalHashBitsStaysInRange() {
		for (int bitsPerKey = 1; bitsPerKey <= 64; bitsPerKey++) {
			int hashBits = SingleWordBloomFilter.optimalHashBits(bitsPerKey);
			assertTrue("bitsPerKey=" + bitsPerKey, hashBits >= 1 && hashBits <= SingleWordBloomFilter.MAX_HASH_BITS);
		}
	}

	@Test
	public void testHashBitsBeyondMaximumRejected() {
		try {
			new SingleWordBloomFilter(SIZE, 16, SingleWordBloomFilter.MAX_HASH_BITS + 1);
			fail("expected IllegalArgumentException for hashBits beyond the maximum");
		} catch (IllegalArgumentException expected) {
			// Expected: beyond that the 64-bit hash runs out of bit positions.
		}
	}

	@Test
	public void testSerializationRoundTrip() throws IOException, ClassNotFoundException {
		// The bucket mask is transient and rebuilt in readObject(); a large-backed filter would address
		// the wrong words if that were missed.
		SingleWordBloomFilter filter = SingleWordBloomFilter.newLargeBacked(SIZE, 16);
		long[] keys = randomKeys(SIZE, 13);
		for (long key : keys) {
			filter.putLong(key);
		}
		ByteArrayOutputStream bytes = new ByteArrayOutputStream();
		try (ObjectOutputStream out = new ObjectOutputStream(bytes)) {
			out.writeObject(filter);
		}
		SingleWordBloomFilter loaded;
		try (ObjectInputStream in = new ObjectInputStream(new ByteArrayInputStream(bytes.toByteArray()))) {
			loaded = (SingleWordBloomFilter) in.readObject();
		}
		assertEquals(filter.getBitSize(), loaded.getBitSize());
		assertEquals(filter.getHashBits(), loaded.getHashBits());
		for (long key : keys) {
			assertTrue(loaded.containsLong(key));
		}
	}

	@Test
	public void testClear() {
		SingleWordBloomFilter filter = new SingleWordBloomFilter(SIZE, 24);
		long[] keys = randomKeys(1000, 5);
		for (long key : keys) {
			filter.putLong(key);
		}
		filter.clear();
		int present = 0;
		for (long key : keys) {
			if (filter.containsLong(key)) {
				present++;
			}
		}
		assertTrue("cleared filter still holds " + present + " keys", present < keys.length / 10);
	}

	@Test
	public void testConcurrentPutMatchesSequentialPut() throws InterruptedException {
		// putLong() locks the backing of the word it touches, so a concurrent fill must end up with
		// exactly the same words as a single-threaded one - ORing is order-independent.
		long[] keys = randomKeys(SIZE, 21);
		for (boolean large : new boolean[] { false, true }) {
			SingleWordBloomFilter sequential = large ? SingleWordBloomFilter.newLargeBacked(SIZE, 16)
					: new SingleWordBloomFilter(SIZE, 16);
			for (long key : keys) {
				sequential.putLong(key);
			}

			SingleWordBloomFilter concurrent = large ? SingleWordBloomFilter.newLargeBacked(SIZE, 16)
					: new SingleWordBloomFilter(SIZE, 16);
			int threads = 8;
			CountDownLatch start = new CountDownLatch(1);
			List<Thread> workers = new ArrayList<Thread>();
			for (int t = 0; t < threads; t++) {
				final int offset = t;
				Thread worker = new Thread(() -> {
					try {
						start.await();
					} catch (InterruptedException e) {
						throw new RuntimeException(e);
					}
					for (int i = offset; i < keys.length; i += threads) {
						concurrent.putLong(keys[i]);
					}
				});
				workers.add(worker);
				worker.start();
			}
			start.countDown();
			for (Thread worker : workers) {
				worker.join();
			}

			long words = sequential.getBitSize() / 64;
			for (long index = 0; index < words; index++) {
				assertEquals("large=" + large + " word " + index, sequential.getWord(index),
						concurrent.getWord(index));
			}
		}
	}
}
