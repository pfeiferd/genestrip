package org.metagene.genestrip.bloom;

import java.util.Random;

import org.metagene.genestrip.trie.KMerTrie;
import org.metagene.genestrip.trie.KMerTrie.KMerTrieVisitor;

import junit.framework.TestCase;

public class CGATBloomFilterTest extends TestCase {
	public void testBloomFilter() {
		int k = 35;
		int size = 5 * 1000;
		double fpp = 0.01;
		
		CGATBloomFilter filter = new CGATBloomFilter(k, size, fpp);
		KMerTrie<String> trie = new KMerTrie<String>(k);
		
		byte[] read = new byte[trie.getLen()];
		byte[] cgat = { 'C', 'G', 'A', 'T' };
		Random random = new Random(42);

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
				assertTrue(filter.contains(read, 0, false));
			}
		}, false);

		trie.visit(new KMerTrieVisitor<String>() {
			@Override
			public void nextValue(KMerTrie<String> trie, byte[] kmer, String value) {
				assertTrue(filter.contains(read, 0, true));
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
		}
		assertTrue(((double) err) / size <= fpp);
	}
}
