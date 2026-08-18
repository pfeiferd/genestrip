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
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import org.junit.Test;
import org.metagene.genestrip.store.KMerStore.UpdateValueProvider;
import org.metagene.genestrip.util.CGAT;

/**
 * Runs the shared {@link AbstractKMerStoreTest} suite against {@link RadixKMerStore} and adds the
 * radix-specific tests (the empty-bucket short circuit and the configurable radix width).
 */
public class RadixKMerStoreTest extends AbstractKMerStoreTest {
	private static final byte[] DECODE_TABLE = CGAT.newDecodeTable();
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
		// radixBits must be >= MIN_RADIX_BITS (16); exercise a few valid widths around the default.
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
	public void testUpdateBatchEqualsUpdate() {
		// Two identical stores: one updated k-mer by k-mer, the other via updateBatch with a small
		// capacity so the batch is flushed several times, ending on a partial flush - exactly the
		// pattern the update reader uses (flush when full and at each region boundary). Both must end
		// up with identical contents. A repeated k-mer checks that duplicates within a batch compose.
		Map<Long, Integer> kmerMap = new LinkedHashMap<Long, Integer>();
		generate(SMALL_SIZE, null, kmerMap);
		long[] kmers = kmerArray(kmerMap);

		RadixKMerStore<Integer> seq = (RadixKMerStore<Integer>) createKMerStore(Integer.class, k, kmers);
		fill(seq, kmerMap);
		seq.optimize();
		RadixKMerStore<Integer> bat = (RadixKMerStore<Integer>) createKMerStore(Integer.class, k, kmers);
		fill(bat, kmerMap);
		bat.optimize();

		// Move every stored value to 1 (idempotent: 1 -> 1 is a no-op); values 0 and 1 always exist.
		UpdateValueProvider<Integer> toOne = oldValue -> 1;

		long[] queries = new long[kmers.length + 1];
		System.arraycopy(kmers, 0, queries, 0, kmers.length);
		queries[kmers.length] = kmers[0]; // duplicate

		for (long kmer : queries) {
			seq.update(kmer, toOne);
		}
		RadixKMerStore.BatchBuffers buf = new RadixKMerStore.BatchBuffers(8);
		for (long kmer : queries) {
			if (buf.add(kmer)) {
				bat.updateBatch(buf, toOne);
			}
		}
		if (!buf.isEmpty()) {
			bat.updateBatch(buf, toOne);
		}

		for (long kmer : kmers) {
			assertEquals(seq.getLong(kmer, null), bat.getLong(kmer, null));
		}
	}

	@Test
	public void testGetBatchEqualsGetLong() {
		// The batched lookup must report exactly what a k-mer by k-mer getLong() reports: the same
		// values, for the same k-mers, in insertion order, with the payload each k-mer was added with.
		// Queries mix stored k-mers, absent ones and a duplicate, and several batch capacities are used
		// so that flushes land on full as well as on partial batches.
		Map<Long, Integer> kmerMap = new LinkedHashMap<Long, Integer>();
		generate(SMALL_SIZE, null, kmerMap);
		long[] kmers = kmerArray(kmerMap);

		RadixKMerStore<Integer> store = (RadixKMerStore<Integer>) createKMerStore(Integer.class, k, kmers);
		fill(store, kmerMap);
		store.optimize();

		// Every second query is a k-mer that is (almost certainly) absent, so misses are covered too.
		long[] queries = new long[2 * kmers.length + 1];
		for (int i = 0; i < kmers.length; i++) {
			queries[2 * i] = kmers[i];
			queries[2 * i + 1] = ~kmers[i];
		}
		queries[queries.length - 1] = kmers[0]; // duplicate

		List<String> expected = new ArrayList<String>();
		for (int i = 0; i < queries.length; i++) {
			Integer value = store.getLong(queries[i], null);
			if (value != null) {
				expected.add(queries[i] + "/" + i + "/" + value);
			}
		}
		assertFalse("expected some hits", expected.isEmpty());

		for (int capacity : new int[] { 1, 8, 128, queries.length + 10 }) {
			List<String> actual = new ArrayList<String>();
			RadixKMerStore.BatchValueConsumer<Integer> consumer = (kmer, payload, value) -> actual
					.add(kmer + "/" + payload + "/" + value);
			RadixKMerStore.BatchBuffers buf = new RadixKMerStore.BatchBuffers(capacity);
			for (int i = 0; i < queries.length; i++) {
				// The index travels as the payload, mirroring how the FT index goal carries its leaf node.
				if (buf.add(queries[i], i)) {
					store.getBatch(buf, consumer);
				}
			}
			if (!buf.isEmpty()) {
				store.getBatch(buf, consumer);
			}
			assertEquals("capacity=" + capacity, expected, actual);
		}
	}

	@Test
	public void testRadixBitsBelowMinimum() {
		// radixBits below MIN_RADIX_BITS (16) is rejected.
		int tooSmall = RadixKMerStore.MIN_RADIX_BITS - 1;
		try {
			new RadixKMerStore<Integer>(k, tooSmall, new int[1 << tooSmall], 0.000001, 0.000001, null, true);
			fail("expected IllegalArgumentException for radixBits below the minimum");
		} catch (IllegalArgumentException expected) {
			// ok
		}
	}

	// Exercises value indices beyond a short's range (65535), up against the per-radix MAX_VALUES,
	// to verify the widened value field actually delivers the expanded capacity and round-trips. The
	// entry's top bit is reserved as the visited mark (see RadixKMerStore#setMarkVisited), so value
	// indices stop below it and a stored entry stays non-negative.
	@Test
	public void testValueCapacityBeyondSortedArray() {
		// The value count is derived from the store's actual per-radix capacity (which depends on
		// RADIX_BITS) rather than hard-coded, so the test stays valid if the default radix width - and
		// hence maxValuesForRadix() - changes. Half-way between the top-bit threshold and the cap, it
		// is guaranteed to satisfy all three bounds below.
		int remainingBits = RadixKMerStore.remainingBitsForRadix(RADIX_BITS);
		int maxValues = RadixKMerStore.maxValuesForRadix(RADIX_BITS);
		// The value field must stop below the reserved mark bit, whatever the radix width is.
		assertEquals("the top entry bit must stay reserved for the visited mark", 1 << (63 - remainingBits),
				maxValues);
		int n = maxValues - maxValues / 4;
		assertTrue(n > Short.MAX_VALUE && n < maxValues);

		// n distinct k-mers, each with a distinct value -> value index == insertion order.
		Map<Long, Integer> kmerMap = new LinkedHashMap<Long, Integer>();
		byte[] read = new byte[k];
		int v = 0;
		while (kmerMap.size() < n) {
			for (int j = 0; j < k; j++) {
				read[j] = DECODE_TABLE[random.nextInt(4)];
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
		// putLong resolves values with a lock-free read, so register them up front (idempotent).
		for (Integer value : kmerMap.values()) {
			store.getAddValueIndex(value);
		}
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
