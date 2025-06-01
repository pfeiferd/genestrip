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
import java.util.List;
import java.util.Map;
import java.util.Random;

import org.junit.Test;
import org.metagene.genestrip.store.KMerStore.KMerStoreVisitor;
import org.metagene.genestrip.util.CGAT;
import org.metagene.genestrip.util.CGATRingBuffer;

import junit.framework.TestCase;

public abstract class AbstractKMerStoreTest extends TestCase {
	protected final Random random = new Random(42);

	protected int k = 31;
	protected int testSize = 1 * 1000 * 1000;
	protected int negativeTestSize = testSize;

	protected abstract <V extends Serializable> KMerStore<V> createKMerStore(Class<V> clazz, Object... params);

	@Test
	public void testSaveLoad() {
		KMerStore<Integer> store = createKMerStore(Integer.class, k);
		store.initSize(testSize);

		Map<List<Byte>, Integer> controlMap = new HashMap<List<Byte>, Integer>();
		Map<Long, Integer> kmerMap = new HashMap<Long, Integer>();
		fillStore(testSize, store, controlMap, kmerMap);
		store.optimize();

		try {
			// Check serialization:
			File tmpdir = Files.createTempDirectory(Paths.get("target"), "serialization-test").toFile();
			tmpdir.deleteOnExit();
			File saved = new File(tmpdir, "saved_trie.ser");
			saved.deleteOnExit();
			if (saved.exists()) {
				saved.delete();
			}
			KMerStore.save(store, saved);

			KMerStore<Integer> loadedStore = KMerStore.load(saved);
			checkStoreContent(testSize, loadedStore, controlMap);
			if (store instanceof KMerSortedArray) {
				checkVisitation(loadedStore, kmerMap);
			}
		} catch (ClassNotFoundException e) {
			throw new RuntimeException(e);
		} catch (IOException e) {
			throw new RuntimeException(e);
		}

	}

	protected void fillStore(int nEntries, KMerStore<Integer> store, Map<List<Byte>, Integer> controlMap,
			Map<Long, Integer> controlMap2) {
		byte[] read = new byte[store.getK()];
		List<Byte> readAsList = null;
		for (int i = 0; i < nEntries; i++) {
			if (controlMap != null) {
				readAsList = new ArrayList<Byte>();
			}
			for (int j = 0; j < store.getK(); j++) {
				read[j] = CGAT.DECODE_TABLE[random.nextInt(4)];
				if (controlMap != null) {
					readAsList.add(read[j]);
				}
			}
			int v = i % KMerSortedArray.MAX_VALUES;
			store.put(read, 0, v);
			if (controlMap != null) {
				controlMap.put(readAsList, v);
			}
			if (controlMap2 != null) {
				controlMap2.put(CGAT.kMerToLong(read, 0, k, null), v);
			}
		}
	}

	protected void checkStoreContent(int negTestValues, KMerStore<Integer> store, Map<List<Byte>, Integer> controlMap) {
		byte[] read = new byte[store.getK()];
		CGATRingBuffer ringBuffer = new CGATRingBuffer(read.length);

		for (List<Byte> key : controlMap.keySet()) {
			Integer v = controlMap.get(key);
			for (int j = 0; j < read.length; j++) {
				read[j] = key.get(j);
				ringBuffer.put(read[j]);
			}

			assertEquals(v, store.get(read, 0));

			CGAT.reverse(read);
			for (int j = 0; j < read.length; j++) {
				ringBuffer.put(read[j]);
			}
			assertEquals(v, store.get(read, 0));
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
				assertNull(store.get(read, 0));
			}

			CGAT.reverse(read);
			readAsList.clear();
			for (int j = 0; j < read.length; j++) {
				readAsList.add(read[j]);
				ringBuffer.put(read[j]);
			}
			if (!controlMap.containsKey(readAsList)) {
				assertNull(store.get(read, 0));
			}
		}
	}

	@Test
	public void testPutGet() {
		KMerStore<Integer> store = createKMerStore(Integer.class, k);
		store.initSize(testSize);

		Map<List<Byte>, Integer> controlMap = new HashMap<List<Byte>, Integer>();
		fillStore(testSize, store, controlMap, null);
		store.optimize();

		checkStoreContent(testSize, store, controlMap);
	}

	@Test
	public void testVisit() {
		KMerStore<Integer> store = createKMerStore(Integer.class, k);
		store.initSize(testSize);

		Map<Long, Integer> kmerMap = new HashMap<Long, Integer>();
		fillStore(testSize, store, null, kmerMap);
		if (store instanceof KMerSortedArray) {
			store.optimize();
		}

		checkVisitation(store, kmerMap);
	}

	protected void checkVisitation(KMerStore<Integer> store, Map<Long, Integer> kmerMap) {
		store.visit(new KMerStoreVisitor<Integer>() {
			@Override
			public void nextValue(KMerStore<Integer> trie, long kmer, Integer value) {
				assertNotNull(value);
				assertEquals(value, kmerMap.get(kmer));
				kmerMap.remove(kmer);
			}
		});
		assertTrue(kmerMap.isEmpty());
	}
}
