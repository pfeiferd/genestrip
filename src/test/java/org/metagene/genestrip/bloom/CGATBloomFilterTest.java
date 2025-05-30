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

import static org.junit.Assert.assertTrue;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import org.junit.Test;
import org.metagene.genestrip.store.KMerTrie;
import org.metagene.genestrip.store.KMerTrie.KMerTrieVisitor;
import org.metagene.genestrip.util.CGAT;

public class CGATBloomFilterTest {
	protected Random random = new Random(42);
	protected int k = 31;
	protected int size = 5 * 1000 * 100;
	protected double fpp = 0.0001;

	@Test
	public void testPutContains() {
		MurmurCGATBloomFilter filter = createFilter(k, size, fpp);

		byte[] reverseRead = new byte[k];
		List<byte[]> reads = new ArrayList<byte[]>();

		for (int i = 1; i < size; i++) {
			byte[] read = new byte[k];
			for (int j = 0; j < k; j++) {
				read[j] = CGAT.DECODE_TABLE[random.nextInt(4)];
				reverseRead[k - j - 1] = CGAT.toComplement(read[j]);
			}
			reads.add(read);

			filter.putLong(CGAT.kMerToLongStraight(read, 0, k, null));
			assertTrue(filter.contains(read, 0, null));
			assertTrue(filter.contains(reverseRead, 0, null));
		}

		for (byte[] read : reads) {
			for (int j = 0; j < k; j++) {
				reverseRead[k - j - 1] = CGAT.toComplement(read[j]);
			}
			assertTrue(filter.contains(read, 0, null));
			assertTrue(filter.contains(reverseRead, 0, null));
		}

		int err = 0;
		byte[] read = new byte[k];
		for (int i = 1; i < size; i++) {
			for (int j = 0; j < k; j++) {
				read[j] = CGAT.DECODE_TABLE[random.nextInt(4)];
			}
			if (filter.contains(read, 0, null)) {
				err++;
			}
		}
		System.out.println("Errors: " + err);
		double testedFp = ((double) err) / (2 * size);
		System.out.println("Tested FP: " + (testedFp * 100) + "%");
		assertTrue(testedFp <= fpp * 1.1);
	}

	@Test
	public void testPutGetViaTrie() {
		MurmurCGATBloomFilter filter = createFilter(k, size, fpp);
		KMerTrie<Integer> trie = new KMerTrie<Integer>(2, k, false);

		byte[] read = new byte[trie.getK()];

		for (int i = 1; i < size; i++) {
			for (int j = 0; j < k; j++) {
				read[j] = CGAT.DECODE_TABLE[random.nextInt(4)];
			}
			filter.putLong(CGAT.kMerToLongStraight(read, 0, k, null));
			trie.put(read, 0, i, false);
		}

		trie.visit(new KMerTrieVisitor<Integer>() {
			@Override
			public void nextValue(KMerTrie<Integer> trie, byte[] kmer, Integer value) {
				assertTrue(filter.contains(kmer, 0, null));
			}
		}, false);

		trie.visit(new KMerTrieVisitor<Integer>() {
			@Override
			public void nextValue(KMerTrie<Integer> trie, byte[] kmer, Integer value) {
				assertTrue(filter.contains(kmer, 0, null));
			}
		}, true);

		int err = 0;
		for (int i = 1; i < size; i++) {
			for (int j = 0; j < k; j++) {
				read[j] = CGAT.DECODE_TABLE[random.nextInt(4)];
			}
			if (filter.contains(read, 0, null)) {
				err++;
			}
		}
		System.out.println("Errors: " + err);
		double testedFp = ((double) err) / (2 * size);
		System.out.println("Tested FP: " + (testedFp * 100) + "%");
		assertTrue(testedFp <= fpp * 1.1);
	}


	protected MurmurCGATBloomFilter createFilter(int k, long size, double fpp) {
		MurmurCGATBloomFilter res = new MurmurCGATBloomFilter(k, fpp);
		res.clear();
		res.ensureExpectedSize(size, isTestLarge());
		return res;
	}
	
	protected boolean isTestLarge() {
		return false;
	}
}
