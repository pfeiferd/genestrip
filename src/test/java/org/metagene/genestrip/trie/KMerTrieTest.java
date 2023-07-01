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

import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;

import org.metagene.genestrip.trie.KMerTrie.KMerTrieVisitor;

public class KMerTrieTest extends AbstractKMerStoreTest {
	public void testVisit2() {
		byte[] cgat = { 'C', 'G', 'A', 'T' };
		Random random = new Random(42);

		for (int k = 0; k < 2; k++) {
			KMerTrie<Integer> trie = createKMerStore(Integer.class, 35);
			Map<List<Byte>, Integer> checkMap = new HashMap<List<Byte>, Integer>();
			byte[] read = new byte[trie.getLen()];
			for (int i = 1; i < 1000; i++) {
				for (int j = 0; j < trie.getLen(); j++) {
					read[j] = cgat[random.nextInt(4)];
				}
				trie.put(read, 0, i, k == 0);

				checkMap.put(byteArrayToList(read), i);
			}

			trie.visit(new KMerTrieVisitor<Integer>() {
				@Override
				public void nextValue(KMerTrie<Integer> trie, byte[] kmer, Integer value) {
					List<Byte> key = byteArrayToList(kmer);
					assertNotNull(value);
					assertEquals(value, checkMap.get(key));
					checkMap.remove(key);
				}
			}, k == 0);
			assertTrue(checkMap.isEmpty());
		}
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
