package org.metagene.genestrip.bloom;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;

import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.Main;
import org.metagene.genestrip.kraken.KrakenResultFastqMergeListener;
import org.metagene.genestrip.kraken.KrakenResultFastqMerger;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;
import org.metagene.genestrip.trie.KMerTrie;
import org.metagene.genestrip.trie.KMerTrie.KMerTrieVisitor;
import org.metagene.genestrip.util.CGAT;
import org.metagene.genestrip.util.CountingDigitTrie;
import org.metagene.genestrip.util.StreamProvider;

import com.google.common.hash.BloomFilter;
import com.google.common.hash.Funnel;
import com.google.common.hash.PrimitiveSink;

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
		int size = 5 * 1000 * 1000;
		double fpp = 0.01;
		
		CGATBloomFilter filter = new CGATBloomFilter(k, size, fpp);
		KMerTrie<Integer> trie = new KMerTrie<Integer>(k);
		
		byte[] read = new byte[trie.getLen()];

		for (int i = 1; i < size; i++) {
			for (int j = 0; j < k; j++) {
				read[j] = cgat[random.nextInt(4)];
			}
			filter.put(read, 0);
			trie.put(read, 0, i, false);
		}
				
		trie.visit(new KMerTrieVisitor<Integer>() {
			@Override
			public void nextValue(KMerTrie<Integer> trie, byte[] kmer, Integer value) {
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
	
	public void testBloomFilterViaProject() throws IOException {
		Main main = new Main();
		main.parseAndRun(new String[] { "bart_h" });

		GSProject project = main.getProject();

		File bartHReads = project.getKmerFastqFile();

		File fromKraken = project.getKrakenOutFile();

		@SuppressWarnings("serial")
		Funnel<byte[]> funnel = new Funnel<byte[]>() {
			@Override
			public void funnel(byte[] from, PrimitiveSink into) {
				for (int i = 0; i < project.getkMserSize(); i++) {
					into.putByte((byte) from[i]);
				}
			}
		};
		
		long size = 5 * 1000 * 1000;
		double fpp = 0.00001;
		
		BloomFilter<byte[]> bloomFilter = BloomFilter.create(funnel, size, fpp);
		CGATBloomFilter cgatBloomFilter = new CGATBloomFilter(project.getkMserSize(), size, fpp);

		@SuppressWarnings("unchecked")
		ObjectGoal<Set<TaxIdNode>, GSProject> taxNodesGoal = (ObjectGoal<Set<TaxIdNode>, GSProject>) main.getGenerator()
				.getGoal("taxids");
		Set<TaxIdNode> nodes = taxNodesGoal.get();

		KrakenResultFastqMergeListener filter = KrakenResultFastqMergeListener.createFilterByTaxIdNodes(nodes,
				new KrakenResultFastqMergeListener() {
					@Override
					public void newTaxIdForRead(long lineCount, byte[] readDescriptor, byte[] read, byte[] readProbs,
							String krakenTaxid, int bps, int pos, String kmerTaxid, int hitLength, byte[] output) {
						cgatBloomFilter.put(read, 0);
						bloomFilter.put(read);
					}
				});

		KrakenResultFastqMerger krakenFilter = new KrakenResultFastqMerger(project.getConfig().getMaxReadSizeBytes());

		InputStream stream1 = StreamProvider.getInputStreamForFile(fromKraken);
		InputStream stream2 = StreamProvider.getInputStreamForFile(bartHReads);
		CountingDigitTrie.print(krakenFilter.process(stream1, stream2, filter), System.out);
		stream1.close();
		stream2.close();

		// Test uncompressed:
		checkFilter(fromKraken, bartHReads, krakenFilter, cgatBloomFilter, bloomFilter, nodes);
	}

	private void checkFilter(File fromKraken, File reads, KrakenResultFastqMerger krakenFilter, CGATBloomFilter filterUnderTest,
			BloomFilter<byte[]> bloomFilter, Set<TaxIdNode> nodes) throws IOException {
		// Positive Test:
		InputStream stream1 = StreamProvider.getInputStreamForFile(fromKraken);
		InputStream stream2 = StreamProvider.getInputStreamForFile(reads);
		krakenFilter.process(stream1, stream2,
				KrakenResultFastqMergeListener.createFilterByTaxIdNodes(nodes, new KrakenResultFastqMergeListener() {
					@Override
					public void newTaxIdForRead(long lineCount, byte[] readDescriptor, byte[] read, byte[] readProbs,
							String krakenTaxid, int bps, int pos, String kmerTaxid, int hitLength, byte[] output) {
						assertTrue(filterUnderTest.contains(read, 0, false));
						assertTrue(bloomFilter.mightContain(read));
					}
				}));
		stream1.close();
		stream2.close();

		// Negative Test:
		byte[] read = new byte[filterUnderTest.getK()];
		byte[] cgat = { 'C', 'G', 'A', 'T' };
		Random random = new Random(42);

		int err = 0;
		int guavaErr = 0;
		for (int i = 1; i < filterUnderTest.getExpectedInsertions(); i++) {
			for (int j = 0; j < read.length; j++) {
				read[j] = cgat[random.nextInt(4)];
			}
			boolean res = filterUnderTest.contains(read, 0, false);
			if (res) {
				err++;
			}
			if (!bloomFilter.mightContain(read)) {
				assertFalse(res);
			}
			else {
				guavaErr++;
			}
		}
		System.out.println("bart_h Errors: " + err);
		double testedFp = ((double) err) / (2 * filterUnderTest.getExpectedInsertions());
		System.out.println("bart_h Tested FP: " + testedFp);
		
		System.out.println("Guava Errors: " + guavaErr);
		double guavaTestedFp = ((double) guavaErr) / (2 * filterUnderTest.getExpectedInsertions());
		System.out.println("bart_h Guavae Tested FP: " + guavaTestedFp);
		
		assertTrue(testedFp <= filterUnderTest.getFpp());
	}
}