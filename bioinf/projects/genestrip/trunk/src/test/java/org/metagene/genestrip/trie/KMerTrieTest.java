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
package org.metagene.genestrip.trie;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.junit.Test;
import org.metagene.genestrip.trie.KMerTrie.KMerTrieVisitor;
import org.metagene.genestrip.util.CGAT;

public class KMerTrieTest extends AbstractKMerStoreTest {
	@Test
	public void testVisitTrie() {
		KMerTrie<Integer> trie = createKMerStore(Integer.class, k);
		trie.initSize(testSize);

		Map<List<Byte>, Integer> checkMap = new HashMap<List<Byte>, Integer>();
		fillStore(testSize, trie, checkMap, null);
		Map<List<Byte>, Integer> checkMap2 = new HashMap<List<Byte>, Integer>(checkMap);

		
		trie.visit(new KMerTrieVisitor<Integer>() {
			@Override
			public void nextValue(KMerTrie<Integer> trie, byte[] kmer, Integer value) {
				List<Byte> key = byteArrayToList(kmer);
				assertNotNull(value);
				assertEquals(value, checkMap.get(key));
				checkMap.remove(key);
			}
		}, false);
		assertTrue(checkMap.isEmpty());

		byte[] copy = new byte[k];
		trie.visit(new KMerTrieVisitor<Integer>() {
			@Override
			public void nextValue(KMerTrie<Integer> trie, byte[] kmer, Integer value) {
				System.arraycopy(kmer, 0, copy, 0, k);
				CGAT.reverse(copy);
				List<Byte> key = byteArrayToList(copy);
				assertNotNull(value);
				assertEquals(value, checkMap2.get(key));
				checkMap2.remove(key);
			}
		}, true);
		assertTrue(checkMap2.isEmpty());
	}

	private List<Byte> byteArrayToList(byte[] arr) {
		List<Byte> res = new ArrayList<Byte>();
		for (int i = 0; i < arr.length; i++) {
			res.add(arr[i]);
		}
		return res;
	}

	@Override
	public <V extends Serializable> KMerTrie<V> createKMerStore(Class<V> clazz, Object... params) {
		return new KMerTrie<V>(2, (int) params[0], false);
	}
}
