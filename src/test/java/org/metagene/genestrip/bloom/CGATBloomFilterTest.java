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

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.Set;

import org.metagene.genestrip.kraken.KrakenResultFastqMergeListener;
import org.metagene.genestrip.kraken.KrakenResultFastqMerger;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;
import org.metagene.genestrip.trie.KMerTrie;
import org.metagene.genestrip.trie.KMerTrie.KMerTrieVisitor;
import org.metagene.genestrip.util.CGAT;
import org.metagene.genestrip.util.StreamProvider;

import junit.framework.TestCase;

public class CGATBloomFilterTest extends TestCase {
	private byte[] cgat = { 'C', 'G', 'A', 'T' };
	private Random random = new Random(42);

	public void testBloomFilter2() {
		int k = 35;
		int size = 5 * 1000 * 1000;
		double fpp = 0.001;

		AbstractCGATBloomFilter filter = createFilter(k, size, fpp);

		byte[] reverseRead = new byte[k];
		List<byte[]> reads = new ArrayList<byte[]>();

		for (int i = 1; i < size; i++) {
			byte[] read = new byte[k];
			for (int j = 0; j < k; j++) {
				read[j] = cgat[random.nextInt(4)];
				reverseRead[k - j - 1] = CGAT.toComplement(read[j]);
			}
			reads.add(read);

			filter.put(read, 0, null);
			assertTrue(filter.contains(read, 0, false, null));
			assertTrue(filter.contains(reverseRead, 0, true, null));
		}

		for (byte[] read : reads) {
			for (int j = 0; j < k; j++) {
				reverseRead[k - j - 1] = CGAT.toComplement(read[j]);
			}
			assertTrue(filter.contains(read, 0, false, null));
			assertTrue(filter.contains(reverseRead, 0, true, null));
		}

		int err = 0;
		byte[] read = new byte[k];
		for (int i = 1; i < size; i++) {
			for (int j = 0; j < k; j++) {
				read[j] = cgat[random.nextInt(4)];
			}
			if (filter.contains(read, 0, false, null)) {
				err++;
			}
			if (filter.contains(read, 0, true, null)) {
				err++;
			}
		}
		System.out.println("Errors: " + err);
		double testedFp = ((double) err) / (2 * size);
		System.out.println("Tested FP: " + (testedFp * 100) + "%");
		assertTrue(testedFp <= fpp * 1.1);
	}

	public void testBloomFilter() {
		int k = 35;
		int size = 5 * 100000;
		double fpp = 0.001;

		AbstractCGATBloomFilter filter = createFilter(k, size, fpp);
		KMerTrie<Integer> trie = new KMerTrie<Integer>(k);

		byte[] read = new byte[trie.getLen()];

		for (int i = 1; i < size; i++) {
			for (int j = 0; j < k; j++) {
				read[j] = cgat[random.nextInt(4)];
			}
			filter.put(read, 0, null);
			trie.put(read, 0, i, false);
		}

		trie.visit(new KMerTrieVisitor<Integer>() {
			@Override
			public void nextValue(KMerTrie<Integer> trie, byte[] kmer, Integer value) {
				assertTrue(filter.contains(kmer, 0, false, null));
			}
		}, false);

		trie.visit(new KMerTrieVisitor<Integer>() {
			@Override
			public void nextValue(KMerTrie<Integer> trie, byte[] kmer, Integer value) {
				assertTrue(filter.contains(kmer, 0, true, null));
			}
		}, true);

		int err = 0;
		for (int i = 1; i < size; i++) {
			for (int j = 0; j < k; j++) {
				read[j] = cgat[random.nextInt(4)];
			}
			if (filter.contains(read, 0, false, null)) {
				err++;
			}
			if (filter.contains(read, 0, true, null)) {
				err++;
			}
		}
		System.out.println("Errors: " + err);
		double testedFp = ((double) err) / (2 * size);
		System.out.println("Tested FP: " + (testedFp * 100) + "%");
		assertTrue(testedFp <= fpp * 1.1);
	}

//	public void testBloomFilterViaProject() throws IOException {
//		Main main = new Main();
//		main.parseAndRun(new String[] { "bart_h", "clear", "genall" });
//
//		GSProject project = main.getProject();
//
//		File bartHReads = ((FileGoal<GSProject>) main.getMaker().getGoal("kmerfastq")).getOutputFile();		
//		File fromKraken = ((FileGoal<GSProject>) main.getMaker().getGoal("kmerkrakenout")).getOutputFile();
//
//		long size = 5 * 1000 * 1000;
//		double fpp = 0.00001;
//
//		AbstractCGATBloomFilter cgatBloomFilter = createFilter(project.getkMserSize(), size, fpp);
//
//		@SuppressWarnings("unchecked")
//		ObjectGoal<Set<TaxIdNode>, GSProject> taxNodesGoal = (ObjectGoal<Set<TaxIdNode>, GSProject>) main.getMaker()
//				.getGoal("taxids");
//		Set<TaxIdNode> nodes = taxNodesGoal.get();
//
//		KrakenResultFastqMergeListener filter = KrakenResultFastqMergeListener.createFilterByTaxIdNodes(nodes,
//				new KrakenResultFastqMergeListener() {
//					@Override
//					public void newTaxIdForRead(long lineCount, byte[] readDescriptor, byte[] read, byte[] readProbs,
//							String krakenTaxid, int bps, int pos, String kmerTaxid, int hitLength, byte[] output) {
//						cgatBloomFilter.put(read, 0, null);
//					}
//				});
//
//		KrakenResultFastqMerger krakenFilter = new KrakenResultFastqMerger(project.getConfig().getMaxReadSizeBytes());
//
//		InputStream stream1 = StreamProvider.getInputStreamForFile(fromKraken);
//		InputStream stream2 = StreamProvider.getInputStreamForFile(bartHReads);
//		CountingDigitTrie.print(krakenFilter.process(stream1, stream2, filter), System.out);
//		stream1.close();
//		stream2.close();
//
//		// Test uncompressed:
//		checkFilter(fromKraken, bartHReads, krakenFilter, cgatBloomFilter, nodes);
//	}

	private void checkFilter(File fromKraken, File reads, KrakenResultFastqMerger krakenFilter,
			AbstractCGATBloomFilter filterUnderTest, Set<TaxIdNode> nodes) throws IOException {
		// Positive Test:
		InputStream stream1 = StreamProvider.getInputStreamForFile(fromKraken);
		InputStream stream2 = StreamProvider.getInputStreamForFile(reads);
		krakenFilter.process(stream1, stream2,
				KrakenResultFastqMergeListener.createFilterByTaxIdNodes(nodes, new KrakenResultFastqMergeListener() {
					@Override
					public void newTaxIdForRead(long lineCount, byte[] readDescriptor, byte[] read, byte[] readProbs,
							String krakenTaxid, int bps, int pos, String kmerTaxid, int hitLength, byte[] output) {
						assertTrue(filterUnderTest.contains(read, 0, false, null));
					}
				}));
		stream1.close();
		stream2.close();

		// Negative Test:
		byte[] read = new byte[filterUnderTest.getK()];
		byte[] cgat = { 'C', 'G', 'A', 'T' };
		Random random = new Random(42);

		int err = 0;
		for (int i = 1; i < filterUnderTest.getExpectedInsertions(); i++) {
			for (int j = 0; j < read.length; j++) {
				read[j] = cgat[random.nextInt(4)];
			}
			if (filterUnderTest.contains(read, 0, false, null)) {
				err++;
			}
		}
		System.out.println("bart_h Errors: " + err);
		double testedFp = ((double) err) / (2 * filterUnderTest.getExpectedInsertions());
		System.out.println("bart_h Tested FP: " + testedFp);

		assertTrue(testedFp <= filterUnderTest.getFpp());
	}
	
	protected AbstractCGATBloomFilter createFilter(int k, long size, double fpp) {
		return new MurmurCGATBloomFilter(k, size, fpp);
	}
}
