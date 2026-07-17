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

import java.io.File;
import java.io.IOException;
import java.io.Serializable;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;

import org.junit.Test;
import org.metagene.genestrip.store.KMerStore.IndexedKMerStoreVisitor;
import org.metagene.genestrip.store.KMerStore.UpdateValueProvider;
import org.metagene.genestrip.store.KMerStore.ValueConverter;
import org.metagene.genestrip.util.CGAT;
import org.metagene.genestrip.util.CGATRingBuffer;

import junit.framework.TestCase;

public abstract class AbstractKMerStoreTest extends TestCase {
	protected final Random random = new Random(42);

	protected int k = 31;
	protected int testSize = 1 * 1000 * 1000;
	protected int negativeTestSize = testSize;

	/**
	 * Creates an empty store for k-mer size {@code k}, sized to hold exactly the given distinct
	 * (canonical) k-mers. Implementations that reserve capacity per radix bucket (e.g.
	 * {@link RadixKMerStore}) derive the sizing from {@code kmers}; the others only need
	 * {@code kmers.length}.
	 */
	protected abstract <V extends Serializable> KMerStore<V> createKMerStore(Class<V> clazz, int k, long[] kmers);

	// --- Shared data generation / store building ------------------------------

	// Generates nEntries random k-mers with values. When non-null, controlMap is filled with
	// (read bytes -> value) (used for the reverse-complement lookups in checkStoreContent) and
	// kmerMap with (canonical k-mer -> value) (used to size + fill the store and to verify
	// visitation). Deterministic via the seeded 'random'.
	protected void generate(int nEntries, Map<List<Byte>, Integer> controlMap, Map<Long, Integer> kmerMap) {
		byte[] read = new byte[k];
		for (int i = 0; i < nEntries; i++) {
			List<Byte> readAsList = controlMap != null ? new ArrayList<Byte>() : null;
			for (int j = 0; j < k; j++) {
				read[j] = CGAT.DECODE_TABLE[random.nextInt(4)];
				if (readAsList != null) {
					readAsList.add(read[j]);
				}
			}
			int v = i % KMerSortedArray.MAX_VALUES;
			if (controlMap != null) {
				controlMap.put(readAsList, v);
			}
			if (kmerMap != null) {
				kmerMap.put(CGAT.kMerToLong(read, 0, k, null), v);
			}
		}
	}

	protected void fill(KMerStore<Integer> store, Map<Long, Integer> kmerMap) {
		// putLong resolves values with a lock-free read, so every value must be registered up front
		// (as the real fill does via the store's initialValues). getAddValueIndex is idempotent.
		for (Integer value : kmerMap.values()) {
			store.getAddValueIndex(value);
		}
		for (Map.Entry<Long, Integer> e : kmerMap.entrySet()) {
			store.putLong(e.getKey(), e.getValue());
		}
	}

	protected static long[] kmerArray(Map<Long, ?> kmerMap) {
		long[] kmers = new long[kmerMap.size()];
		int i = 0;
		for (long kmer : kmerMap.keySet()) {
			kmers[i++] = kmer;
		}
		return kmers;
	}

	/** Generates {@code nEntries} k-mers, creates a store sized for them, fills and optimizes it. */
	protected KMerStore<Integer> buildStore(int nEntries, Map<List<Byte>, Integer> controlMap,
			Map<Long, Integer> kmerMap) {
		generate(nEntries, controlMap, kmerMap);
		KMerStore<Integer> store = createKMerStore(Integer.class, k, kmerArray(kmerMap));
		fill(store, kmerMap);
		store.optimize();
		return store;
	}

	// --- Shared tests ---------------------------------------------------------

	@Test
	public void testPutGet() {
		Map<List<Byte>, Integer> controlMap = new HashMap<List<Byte>, Integer>();
		KMerStore<Integer> store = buildStore(testSize, controlMap, new LinkedHashMap<Long, Integer>());
		checkStoreContent(store, controlMap);
	}

	@Test
	public void testVisit() {
		Map<Long, Integer> kmerMap = new LinkedHashMap<Long, Integer>();
		KMerStore<Integer> store = buildStore(testSize, null, kmerMap);
		checkVisitation(store, kmerMap);
	}

	@Test
	public void testSaveLoad() {
		Map<List<Byte>, Integer> controlMap = new HashMap<List<Byte>, Integer>();
		Map<Long, Integer> kmerMap = new LinkedHashMap<Long, Integer>();
		KMerStore<Integer> store = buildStore(testSize, controlMap, kmerMap);

		try {
			File tmpdir = Files.createTempDirectory(Paths.get("target"), "serialization-test").toFile();
			tmpdir.deleteOnExit();
			File saved = new File(tmpdir, "saved_store.ser");
			saved.deleteOnExit();
			if (saved.exists()) {
				saved.delete();
			}
			KMerStore.save(store, saved);

			KMerStore<Integer> loadedStore = KMerStore.load(saved);
			checkStoreContent(loadedStore, controlMap);
			checkVisitation(loadedStore, kmerMap);
		} catch (ClassNotFoundException e) {
			throw new RuntimeException(e);
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}

	@Test
	public void testUpdate() {
		Map<Long, Integer> kmerMap = new LinkedHashMap<Long, Integer>();
		KMerStore<Integer> store = buildStore(testSize, null, kmerMap);

		Map.Entry<Long, Integer> first = kmerMap.entrySet().iterator().next();
		long kmer = first.getKey();
		// Update to a different value that is already present (values 0 and 1 always are), so no new
		// value index is created even when the store is filled up to its MAX_VALUES distinct values
		// (as KMerSortedArray is here).
		final int newValue = first.getValue() == 0 ? 1 : 0;
		assertTrue(store.update(kmer, new UpdateValueProvider<Integer>() {
			@Override
			public Integer getUpdateValue(Integer oldValue) {
				return newValue;
			}
		}));
		assertEquals(Integer.valueOf(newValue), store.getLong(kmer, null));
	}

	@Test
	public void testConvertValues() {
		Map<Long, Integer> kmerMap = new LinkedHashMap<Long, Integer>();
		KMerStore<Integer> store = buildStore(testSize, null, kmerMap);

		KMerStore<String> converted = store.convertValues(new ValueConverter<Integer, String>() {
			@Override
			public String convertValue(Integer value) {
				return "v" + value;
			}
		});
		for (Map.Entry<Long, Integer> e : kmerMap.entrySet()) {
			assertEquals("v" + e.getValue(), converted.getLong(e.getKey(), null));
		}
	}

	// --- Shared verification helpers ------------------------------------------

	protected void checkStoreContent(KMerStore<Integer> store, Map<List<Byte>, Integer> controlMap) {
		byte[] read = new byte[store.getK()];
		CGATRingBuffer ringBuffer = new CGATRingBuffer(read.length);
		long[] pos = new long[1];
		long entries = store.getEntries();

		for (List<Byte> key : controlMap.keySet()) {
			Integer v = controlMap.get(key);
			for (int j = 0; j < read.length; j++) {
				read[j] = key.get(j);
				ringBuffer.put(read[j]);
			}

			assertEquals(v, store.getLong(CGAT.kMerToLong(read, 0, store.getK(), null), pos));
			assertTrue("hit position must be within [0, entries)", pos[0] >= 0 && pos[0] < entries);

			CGAT.reverse(read);
			for (int j = 0; j < read.length; j++) {
				ringBuffer.put(read[j]);
			}
			assertEquals(v, store.getLong(CGAT.kMerToLong(read, 0, store.getK(), null), null));
		}

		List<Byte> readAsList = new ArrayList<Byte>();
		for (int i = 1; i <= negativeTestSize; i++) {
			readAsList.clear();
			for (int j = 0; j < read.length; j++) {
				read[j] = CGAT.DECODE_TABLE[random.nextInt(4)];
				readAsList.add(read[j]);
				ringBuffer.put(read[j]);
			}
			if (!controlMap.containsKey(readAsList)) {
				assertNull(store.getLong(CGAT.kMerToLong(read, 0, store.getK(), null), null));
			}

			CGAT.reverse(read);
			readAsList.clear();
			for (int j = 0; j < read.length; j++) {
				readAsList.add(read[j]);
				ringBuffer.put(read[j]);
			}
			if (!controlMap.containsKey(readAsList)) {
				assertNull(store.getLong(CGAT.kMerToLong(read, 0, store.getK(), null), null));
			}
		}
	}

	protected void checkVisitation(KMerStore<Integer> store, Map<Long, Integer> kmerMap) {
		Map<Long, Integer> remaining = new HashMap<Long, Integer>(kmerMap);
		Set<Long> positions = new HashSet<Long>();
		long entries = store.getEntries();
		store.visit(new IndexedKMerStoreVisitor<Integer>() {
			@Override
			public void nextValue(KMerStore<Integer> trie, long kmer, int index, long pos) {
				Integer value = trie.getValueForIndex(index);
				assertNotNull(value);
				assertEquals(value, remaining.remove(kmer));
				assertTrue("visit position must be unique and within [0, entries)",
						pos >= 0 && pos < entries && positions.add(pos));
			}
		});
		assertTrue(remaining.isEmpty());
		assertEquals(entries, (long) positions.size());
	}
}
