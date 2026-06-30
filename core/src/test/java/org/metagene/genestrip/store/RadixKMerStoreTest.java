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

import java.io.Serializable;
import java.util.LinkedHashMap;
import java.util.Map;

import org.junit.Test;
import org.metagene.genestrip.util.CGAT;

/**
 * Runs the shared {@link AbstractKMerStoreTest} suite against {@link RadixKMerStore} and adds the
 * radix-specific tests (the empty-bucket short circuit and the configurable radix width).
 */
public class RadixKMerStoreTest extends AbstractKMerStoreTest {
	private static final int RADIX_BITS = RadixKMerStore.DEFAULT_RADIX_BITS;
	// A k-mer set small enough that many radix buckets stay empty (for the early-out test).
	private static final int SMALL_SIZE = 20000;

	@Override
	public <V extends Serializable> KMerStore<V> createKMerStore(Class<V> clazz, int k, long[] kmers) {
		int[] bucketSizes = new int[1 << RADIX_BITS];
		for (long kmer : kmers) {
			bucketSizes[RadixKMerStore.radixOf(kmer, RADIX_BITS)]++;
		}
		return new RadixKMerStore<V>(k, RADIX_BITS, bucketSizes, 0.000001, 0.000001, null, true);
	}

	@Test
	public void testRadixEarlyOut() {
		Map<Long, Integer> kmerMap = new LinkedHashMap<Long, Integer>();
		KMerStore<Integer> store = buildStore(SMALL_SIZE, null, kmerMap);

		boolean[] used = new boolean[1 << RADIX_BITS];
		for (long kmer : kmerMap.keySet()) {
			used[RadixKMerStore.radixOf(kmer, RADIX_BITS)] = true;
		}
		// A k-mer whose low RADIX_BITS bits address an empty bucket must short-circuit to null.
		// (kmer == radix has zero remaining bits and radixOf(kmer, RADIX_BITS) == radix.)
		int tested = 0;
		for (int radix = 0; radix < used.length && tested < 100; radix++) {
			if (!used[radix]) {
				assertNull("empty radix bucket must return null", store.getLong(radix, null));
				tested++;
			}
		}
		assertTrue("expected some empty radix buckets", tested > 0);
	}

	@Test
	public void testConfigurableRadixBits() {
		Map<Long, Integer> kmerMap = new LinkedHashMap<Long, Integer>();
		generate(SMALL_SIZE, null, kmerMap);
		long[] kmers = kmerArray(kmerMap);
		// radixBits must be >= MIN_RADIX_BITS (17); exercise a few valid widths around the default.
		for (int radixBits : new int[] { RadixKMerStore.MIN_RADIX_BITS, 20, 22 }) {
			int[] sizes = new int[1 << radixBits];
			for (long kmer : kmers) {
				sizes[RadixKMerStore.radixOf(kmer, radixBits)]++;
			}
			RadixKMerStore<Integer> store = new RadixKMerStore<Integer>(k, radixBits, sizes, 0.000001, 0.000001, null, true);
			fill(store, kmerMap);
			store.optimize();
			assertEquals(radixBits, store.getRadixBits());
			for (Map.Entry<Long, Integer> e : kmerMap.entrySet()) {
				assertEquals("radixBits=" + radixBits, e.getValue(), store.getLong(e.getKey(), null));
			}
		}
	}

	@Test
	public void testRadixBitsBelowMinimum() {
		// radixBits below MIN_RADIX_BITS (17) is rejected.
		int tooSmall = RadixKMerStore.MIN_RADIX_BITS - 1;
		try {
			new RadixKMerStore<Integer>(k, tooSmall, new int[1 << tooSmall], 0.000001, 0.000001, null, true);
			fail("expected IllegalArgumentException for radixBits below the minimum");
		} catch (IllegalArgumentException expected) {
			// ok
		}
	}

	// Exercises value indices beyond KMerSortedArray's cap (65535) and beyond 2^(63-REMAINING_BITS) -
	// i.e. into the range where the packed entry's high bit is set (a negative long) - to verify the
	// 19-bit value field actually delivers the expanded RadixKMerStore.MAX_VALUES and round-trips.
	@Test
	public void testValueCapacityBeyondSortedArray() {
		int n = 300_000;
		assertTrue(n > KMerSortedArray.MAX_VALUES
				&& n > (1 << (63 - RadixKMerStore.REMAINING_BITS)) // some indices set the entry's top bit
				&& n < RadixKMerStore.MAX_VALUES);

		// n distinct k-mers, each with a distinct value -> value index == insertion order.
		Map<Long, Integer> kmerMap = new LinkedHashMap<Long, Integer>();
		byte[] read = new byte[k];
		int v = 0;
		while (kmerMap.size() < n) {
			for (int j = 0; j < k; j++) {
				read[j] = CGAT.DECODE_TABLE[random.nextInt(4)];
			}
			long kmer = CGAT.kMerToLong(read, 0, k, null);
			if (!kmerMap.containsKey(kmer)) {
				kmerMap.put(kmer, v++);
			}
		}

		long[] kmers = kmerArray(kmerMap);
		int[] bucketSizes = new int[1 << RADIX_BITS];
		for (long kmer : kmers) {
			bucketSizes[RadixKMerStore.radixOf(kmer, RADIX_BITS)]++;
		}
		// A very low fpp so no k-mer is dropped on a fill-time filter false positive.
		RadixKMerStore<Integer> store = new RadixKMerStore<Integer>(k, RADIX_BITS, bucketSizes, 1e-9, 1e-9, null, true);
		for (Map.Entry<Long, Integer> e : kmerMap.entrySet()) {
			assertTrue(store.putLong(e.getKey(), e.getValue()));
		}
		store.optimize();

		assertEquals(n, store.getNValues());
		for (Map.Entry<Long, Integer> e : kmerMap.entrySet()) {
			assertEquals(e.getValue(), store.getLong(e.getKey(), null));
		}
	}
}
