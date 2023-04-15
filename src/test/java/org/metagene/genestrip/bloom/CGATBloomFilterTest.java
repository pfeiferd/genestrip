package org.metagene.genestrip.bloom;

import java.util.Random;

import org.metagene.genestrip.trie.KMerTrie;
import org.metagene.genestrip.trie.KMerTrie.KMerTrieVisitor;
import org.metagene.genestrip.util.CGAT;

import junit.framework.TestCase;

public class CGATBloomFilterTest extends TestCase {
	private byte[] cgat = { 'C', 'G', 'A', 'T' };
	private Random random = new Random(42);
	
	public void testHash() {
		int k = 35;
		int size = 5 * 1000;
		double fpp = 0.01;

		CGATBloomFilter filter = new CGATBloomFilter(k, size, fpp);
		
		byte[] read = new byte[k];
		byte[] reverseRead = new byte[k];
		
		for (int i = 1; i < size; i++) {
			for (int j = 0; j < k; j++) {
				read[j] = cgat[random.nextInt(4)];
				reverseRead[k - j - 1] = CGAT.toComplement(read[j]);
			}
			long hash1 = filter.hash(read, i, false, false);
			long hash2 = filter.hash(reverseRead, i, false, true);
			
			assertEquals(hash1, hash2);
			
			filter.put(read, 0);
			assertTrue(filter.contains(reverseRead, 0, true));
		}		
	}
	
	public void testBloomFilter() {
		int k = 35;
		int size = 5 * 1000;
		double fpp = 0.01;
		
		CGATBloomFilter filter = new CGATBloomFilter(k, size, fpp);
		KMerTrie<String> trie = new KMerTrie<String>(k);
		
		byte[] read = new byte[trie.getLen()];

		for (int i = 1; i < size; i++) {
			for (int j = 0; j < k; j++) {
				read[j] = cgat[random.nextInt(4)];
			}
			filter.put(read, 0);
			trie.put(read, 0, "-", false);
		}
		
		trie.visit(new KMerTrieVisitor<String>() {
			@Override
			public void nextValue(KMerTrie<String> trie, byte[] kmer, String value) {
				assertTrue(filter.contains(kmer, 0, false));
			}
		}, false);

		trie.visit(new KMerTrieVisitor<String>() {
			@Override
			public void nextValue(KMerTrie<String> trie, byte[] kmer, String value) {
				assertTrue(filter.contains(kmer, 0, true));
			}
		}, true);
		
		int err = 0;
		for (int i = 1; i < size; i++) {
			for (int j = 0; j < k; j++) {
				read[j] = cgat[random.nextInt(4)];
			}
			if (filter.contains(read, 0, false)) {
				err++;
			}
			if (filter.contains(read, 0, true)) {
				err++;
			}
		}
		System.out.println("Errors: " + err);
		double testedFp = ((double) err) / (2 * size);
		System.out.println("Tested FP: " + testedFp);
		assertTrue(testedFp <= fpp);
	}
}
