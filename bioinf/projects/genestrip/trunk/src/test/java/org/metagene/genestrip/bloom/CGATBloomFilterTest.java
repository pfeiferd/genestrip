package org.metagene.genestrip.bloom;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
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
		int size = 5 * 10000;
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
			assertTrue(filter.contains(read, 0, false));
			assertTrue(filter.contains(reverseRead, 0, true));
		}		
	}
	
	public void testBloomFilter2() {
		int k = 35;
		int size = 5 * 1000 * 1000;
		double fpp = 0.01;

		CGATBloomFilter filter = new CGATBloomFilter(k, size, fpp);
		
		byte[] reverseRead = new byte[k];
		List<byte[]> reads = new ArrayList<byte[]>();
		
		for (int i = 1; i < size; i++) {
			byte[] read = new byte[k];
			for (int j = 0; j < k; j++) {
				read[j] = cgat[random.nextInt(4)];
				reverseRead[k - j - 1] = CGAT.toComplement(read[j]);
			}
			reads.add(read);
			long hash1 = filter.hash(read, i, false, false);
			long hash2 = filter.hash(reverseRead, i, false, true);
			
			assertEquals(hash1, hash2);
			
			filter.put(read, 0);
			assertTrue(filter.contains(read, 0, false));
			assertTrue(filter.contains(reverseRead, 0, true));
		}		
		
		for (byte[] read : reads) {
			for (int j = 0; j < k; j++) {
				reverseRead[k - j - 1] = CGAT.toComplement(read[j]);
			}
			assertTrue(filter.contains(read, 0, false));
			assertTrue(filter.contains(reverseRead, 0, true));
		}		
		
		int err = 0;
		byte[] read = new byte[k];
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
	
	public void testBloomFilter() {
		int k = 35;
		int size = 5 * 1000;
		double fpp = 0.01;
		
		CGATBloomFilter filter = new CGATBloomFilter(k, size, fpp);
		KMerTrie<Integer> trie = new KMerTrie<Integer>(k);
		
		byte[] read = new byte[trie.getLen()];
		Map<List<Byte>, Integer> checkMap = new HashMap<List<Byte>, Integer>();

		for (int i = 1; i < size; i++) {
			for (int j = 0; j < k; j++) {
				read[j] = cgat[random.nextInt(4)];
			}
			filter.put(read, 0);
			trie.put(read, 0, i, false);
			checkMap.put(byteArrayToList(read), i);
		}
		
		trie.visit(new KMerTrieVisitor<Integer>() {
			int count;
			
			@Override
			public void nextValue(KMerTrie<Integer> trie, byte[] kmer, Integer value) {
				for (int i = 0; i < k; i++) {
					assertTrue(CGAT.isCGAT(kmer[i]));
				}
				System.out.println(++count);
				List<Byte> key = byteArrayToList(kmer);
				assertEquals(value, checkMap.get(key));
				assertTrue(filter.contains(kmer, 0, false));
			}
		}, false);

		trie.visit(new KMerTrieVisitor<Integer>() {
			@Override
			public void nextValue(KMerTrie<Integer> trie, byte[] kmer, Integer value) {
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

	private List<Byte> byteArrayToList(byte[] arr) {
		List<Byte> res = new ArrayList<Byte>();
		for (int i = 0; i < arr.length; i++) {
			res.add(arr[i]);
		}
		return res;
	}
}
