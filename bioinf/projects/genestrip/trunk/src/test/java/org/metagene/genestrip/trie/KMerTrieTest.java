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

import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;
import java.util.zip.GZIPInputStream;

import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.Main;
import org.metagene.genestrip.kraken.KrakenResultFastqMergeListener;
import org.metagene.genestrip.kraken.KrakenResultFastqMerger;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;
import org.metagene.genestrip.trie.KMerTrie.KMerTrieVisitor;
import org.metagene.genestrip.util.ByteArrayToString;
import org.metagene.genestrip.util.CountingDigitTrie;

import com.google.common.hash.BloomFilter;
import com.google.common.hash.Funnel;
import com.google.common.hash.PrimitiveSink;

import junit.framework.TestCase;

public class KMerTrieTest extends TestCase {
	public void testPutGet() throws IOException {
		Main main = new Main();
		main.parseAndRun(new String[] { "bart_h" });

		GSProject project = main.getProject();

		File bartHReads = project.getKmerFastqFile();
		FileInputStream fStream = new FileInputStream(bartHReads);
		GZIPInputStream gStream = new GZIPInputStream(fStream, 4096);

		File fromKraken = project.getKrakenOutFile();

		KMerTrie<String> trie = new KMerTrie<String>(project.getkMserSize());

		@SuppressWarnings("serial")
		Funnel<byte[]> funnel = new Funnel<byte[]>() {
			@Override
			public void funnel(byte[] from, PrimitiveSink into) {
				for (int i = 0; i < trie.getLen(); i++) {
					into.putByte((byte) from[i]);
				}
			}
		};
		BloomFilter<byte[]> bloomFilter = BloomFilter.create(funnel, 5 * 1000 * 1000, 0.00001);

		@SuppressWarnings("unchecked")
		ObjectGoal<Set<TaxIdNode>, GSProject> taxNodesGoal = (ObjectGoal<Set<TaxIdNode>, GSProject>) main.getGenerator()
				.getGoal("taxids");
		Set<TaxIdNode> nodes = taxNodesGoal.get();

		KrakenResultFastqMergeListener filter = KrakenResultFastqMergeListener.createFilterByTaxIdNodes(nodes,
				KrakenResultFastqMergeListener.fillKMerTrie(trie, new KrakenResultFastqMergeListener() {
					@Override
					public void newTaxIdForRead(long lineCount, byte[] readDescriptor, byte[] read, byte[] readProbs,
							String krakenTaxid, int bps, int pos, String kmerTaxid, int hitLength, byte[] output) {
						bloomFilter.put(read);
					}
				}));

		KrakenResultFastqMerger krakenFilter = new KrakenResultFastqMerger(project.getConfig().getMaxReadSizeBytes());

		CountingDigitTrie.print(krakenFilter.process(new BufferedInputStream(new FileInputStream(fromKraken)),
				new BufferedInputStream(gStream), filter), System.out);

		// Test uncompressed:
		checkTrie(fromKraken, bartHReads, krakenFilter, trie, bloomFilter, nodes);

		trie.compress();

		// Test compressed;
		checkTrie(fromKraken, bartHReads, krakenFilter, trie, bloomFilter, nodes);
	}

	private void checkTrie(File fromKraken, File reads, KrakenResultFastqMerger krakenFilter, KMerTrie<String> trie,
			BloomFilter<byte[]> bloomFilter, Set<TaxIdNode> nodes) throws IOException {
		FileInputStream fStream2 = new FileInputStream(reads);
		GZIPInputStream gStream2 = new GZIPInputStream(fStream2, 4096);

		// Positive Test:
		krakenFilter.process(new BufferedInputStream(new FileInputStream(fromKraken)),
				new BufferedInputStream(gStream2),
				KrakenResultFastqMergeListener.createFilterByTaxIdNodes(nodes, new KrakenResultFastqMergeListener() {
					@Override
					public void newTaxIdForRead(long lineCount, byte[] readDescriptor, byte[] read, byte[] readProbs,
							String krakenTaxid, int bps, int pos, String kmerTaxid, int hitLength, byte[] output) {
						assertEquals(kmerTaxid, trie.get(read, 0, false));
					}
				}));

		// Negative Test:
		byte[] read = new byte[trie.getLen()];
		byte[] cgat = { 'C', 'G', 'A', 'T' };
		Random random = new Random(42);

		for (int i = 1; i < 5 * 1000 * 1000; i++) {
			for (int j = 0; j < trie.getLen(); j++) {
				read[j] = cgat[random.nextInt(4)];
			}
			if (!bloomFilter.mightContain(read)) {
				assertNull(trie.get(read, 0, false));
			}
		}
	}
	
	public void testVisit() {
		KMerTrie<Integer> trie = new KMerTrie<Integer>(2, 35, false);
		Map<List<Object>, Integer> checkMap = new HashMap<List<Object>, Integer>();

		byte[] cgat = { 'C', 'G', 'A', 'T' };
		Random random = new Random(42);

		for (int i = 1; i < 1000; i++) {
			byte[] read = new byte[trie.getLen()];
			for (int j = 0; j < trie.getLen(); j++) {
				read[j] = cgat[random.nextInt(4)];
			}
			if (i == 356) {
				System.out.println(ByteArrayToString.toString(read));				
			}
			trie.put(read, 0, i);

			checkMap.put(Arrays.asList(read), i);
		}
		
		trie.visit(new KMerTrieVisitor<Integer>() {
			@Override
			public void nextValue(KMerTrie<Integer> trie, byte[] kmer, Integer value) {
				System.out.println(ByteArrayToString.toString(kmer));
				List<Object> key = Arrays.asList(kmer);
				assertEquals(value,checkMap.get(key));
				checkMap.remove(key);
			}
		});
		
		assertTrue(checkMap.isEmpty());
	}
}