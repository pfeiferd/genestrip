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
import org.metagene.genestrip.bloom.AbstractCGATBloomFilter;
import org.metagene.genestrip.bloom.MurmurCGATBloomFilter;
import org.metagene.genestrip.kraken.KrakenResultFastqMergeListener;
import org.metagene.genestrip.kraken.KrakenResultFastqMerger;
import org.metagene.genestrip.make.FileGoal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;
import org.metagene.genestrip.trie.KMerTrie.KMerTrieVisitor;
import org.metagene.genestrip.util.CountingDigitTrie;
import org.metagene.genestrip.util.StreamProvider;

import junit.framework.TestCase;

public class KMerTrieTest extends TestCase {
	public void testPutGet() throws IOException {
		Main main = new Main();
		main.parseAndRun(new String[] { "bart_h", "clear", "genall" });

		GSProject project = main.getProject();

		File bartHReads = ((FileGoal<GSProject>) main.getMaker().getGoal("kmerfastq")).getFile();		
		File fromKraken = ((FileGoal<GSProject>) main.getMaker().getGoal("kmerkrakenout")).getFile();

		KMerTrie<String> trie = new KMerTrie<String>(project.getkMserSize());

		AbstractCGATBloomFilter bloomFilter = new MurmurCGATBloomFilter(trie.getLen(), 5 * 1000 * 1000, 0.00001);

		@SuppressWarnings("unchecked")
		ObjectGoal<Set<TaxIdNode>, GSProject> taxNodesGoal = (ObjectGoal<Set<TaxIdNode>, GSProject>) main.getMaker()
				.getGoal("taxids");
		Set<TaxIdNode> nodes = taxNodesGoal.get();

		KrakenResultFastqMergeListener filter = KrakenResultFastqMergeListener.createFilterByTaxIdNodes(nodes,
				new KrakenResultFastqMergeListener() {
					@Override
					public void newTaxIdForRead(long lineCount, byte[] readDescriptor, byte[] read, byte[] readProbs,
							String krakenTaxid, int bps, int pos, String kmerTaxid, int hitLength, byte[] output) {
						trie.put(read, 0, kmerTaxid, false);
						bloomFilter.put(read, 0, null);
					}
				});

		KrakenResultFastqMerger krakenFilter = new KrakenResultFastqMerger(project.getConfig().getMaxReadSizeBytes());

		InputStream stream1 = StreamProvider.getInputStreamForFile(fromKraken);
		InputStream stream2 = StreamProvider.getInputStreamForFile(bartHReads);
		CountingDigitTrie.print(krakenFilter.process(stream1, stream2, filter), System.out);
		stream1.close();
		stream2.close();

		// Test uncompressed:
		checkTrie(fromKraken, bartHReads, krakenFilter, trie, bloomFilter, nodes);

		trie.compress();

		// Test compressed;
		checkTrie(fromKraken, bartHReads, krakenFilter, trie, bloomFilter, nodes);
	}

	private void checkTrie(File fromKraken, File reads, KrakenResultFastqMerger krakenFilter, KMerTrie<String> trie,
			AbstractCGATBloomFilter bloomFilter, Set<TaxIdNode> nodes) throws IOException {
		// Positive Test:
		InputStream stream1 = StreamProvider.getInputStreamForFile(fromKraken);
		InputStream stream2 = StreamProvider.getInputStreamForFile(reads);
		krakenFilter.process(stream1, stream2,
				KrakenResultFastqMergeListener.createFilterByTaxIdNodes(nodes, new KrakenResultFastqMergeListener() {
					@Override
					public void newTaxIdForRead(long lineCount, byte[] readDescriptor, byte[] read, byte[] readProbs,
							String krakenTaxid, int bps, int pos, String kmerTaxid, int hitLength, byte[] output) {
						assertEquals(kmerTaxid, trie.get(read, 0, false));
					}
				}));
		stream1.close();
		stream2.close();

		// Negative Test:
		byte[] read = new byte[trie.getLen()];
		byte[] cgat = { 'C', 'G', 'A', 'T' };
		Random random = new Random(42);

		for (int i = 1; i < 5 * 1000 * 1000; i++) {
			for (int j = 0; j < trie.getLen(); j++) {
				read[j] = cgat[random.nextInt(4)];
			}
			if (!bloomFilter.contains(read, 0, false, null)) {
				assertNull(trie.get(read, 0, false));
			}
		}
	}

	public void testVisit() {
		byte[] cgat = { 'C', 'G', 'A', 'T' };
		Random random = new Random(42);

		for (int k = 0; k < 2; k++) {
			KMerTrie<Integer> trie = new KMerTrie<Integer>(2, 35, false);
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
}
