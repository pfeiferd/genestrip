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
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.HashMap;
import java.util.Map;
import java.util.Random;
import java.util.Set;

import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.Main;
import org.metagene.genestrip.bloom.MurmurCGATBloomFilter;
import org.metagene.genestrip.goals.KrakenResCountGoal;
import org.metagene.genestrip.kraken.KrakenResultFastqMergeListener;
import org.metagene.genestrip.kraken.KrakenResultFastqMerger;
import org.metagene.genestrip.make.FileGoal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;
import org.metagene.genestrip.trie.KMerStore.KMerStoreVisitor;
import org.metagene.genestrip.util.CGAT;
import org.metagene.genestrip.util.StreamProvider;

import junit.framework.TestCase;

public abstract class AbstractKMerStoreTest extends TestCase implements KMerStoreFactory {
	protected final int testSize = 1000000;
	protected final int negativeTestSize = 5 * 1000 * 1000;
	
	public void testPutGet() throws IOException {
		Main main = new Main();
		main.parseAndRun(new String[] { "bart_h", "clear", "genall" });

		GSProject project = main.getProject();

		File bartHReads = ((FileGoal<GSProject>) main.getMaker().getGoal("kmerfastq")).getFile();

		String outGoal = project.getConfig().isUseKraken1() ? "sort" : "kmerkrakenout";
		File fromKraken = ((FileGoal<GSProject>) main.getMaker().getGoal(outGoal)).getFile();

		final KMerStore<String> trie = createKMerStore(String.class, project.getkMserSize());

		MurmurCGATBloomFilter bloomFilter = new MurmurCGATBloomFilter(trie.getLen(), 0.00001);
		bloomFilter.clearAndEnsureCapacity(5 * 1000 * 1000);

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
						bloomFilter.put(read, 0);
					}
				});

		KrakenResultFastqMerger krakenFilter = new KrakenResultFastqMerger(project.getConfig().getMaxReadSizeBytes());

		InputStream stream1 = StreamProvider.getInputStreamForFile(fromKraken);
		InputStream stream2 = StreamProvider.getInputStreamForFile(bartHReads);
		KrakenResCountGoal.print(krakenFilter.process(stream1, stream2, filter), System.out);
		stream1.close();
		stream2.close();

		// Test uncompressed:
		checkTrie(fromKraken, bartHReads, krakenFilter, trie, bloomFilter, nodes);

		trie.optimize();

		// Test compressed;
		checkTrie(fromKraken, bartHReads, krakenFilter, trie, bloomFilter, nodes);

		// Check serialization:
		File tmpdir = Files.createTempDirectory(Paths.get("target"), "serialization-test").toFile();
		tmpdir.deleteOnExit();
		File saved = new File(tmpdir, "saved_trie.ser");
		saved.deleteOnExit();
		if (saved.exists()) {
			saved.delete();
		}
		KMerStore.save(trie, saved);
		try {
			KMerStore<String> loadedTrie = KMerStore.load(saved);
			checkTrie(fromKraken, bartHReads, krakenFilter, loadedTrie, bloomFilter, nodes);
		} catch (ClassNotFoundException e) {
			throw new RuntimeException(e);
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}

	private void checkTrie(File fromKraken, File reads, KrakenResultFastqMerger krakenFilter, KMerStore<String> trie,
			MurmurCGATBloomFilter bloomFilter, Set<TaxIdNode> nodes) throws IOException {
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

		for (int i = 1; i <= negativeTestSize; i++) {
			for (int j = 0; j < trie.getLen(); j++) {
				read[j] = cgat[random.nextInt(4)];
			}
			if (!bloomFilter.containsStraight(read, 0, null)) {
				assertNull(trie.get(read, 0, false));
			}
		}
	}

	public void testSaveLoad() {
		byte[] cgat = { 'C', 'G', 'A', 'T' };
		Random random = new Random(42);

		KMerStore<Integer> trie = createKMerStore(Integer.class, 35);
		trie.initSize(testSize);
		
		Map<Long, Integer> checkMap = new HashMap<Long, Integer>();
		byte[] read = new byte[trie.getLen()];
		for (int i = 1; i <= testSize; i++) {
			for (int j = 0; j < trie.getLen(); j++) {
				read[j] = cgat[random.nextInt(4)];
			}
			trie.put(read, 0, i, false);
			checkMap.put(CGAT.kmerToLongStraight(read, 0, trie.getLen(), null), i);
		}
		trie.optimize();
		
		try {
			// Check serialization:
			File tmpdir = Files.createTempDirectory(Paths.get("target"), "serialization-test").toFile();
			tmpdir.deleteOnExit();
			File saved = new File(tmpdir, "saved_trie.ser");
			saved.deleteOnExit();
			if (saved.exists()) {
				saved.delete();
			}
			KMerStore.save(trie, saved);

			KMerStore<Integer> loadedTrie = KMerStore.load(saved);
			loadedTrie.visit(new KMerStoreVisitor<Integer>() {
				@Override
				public void nextValue(KMerStore<Integer> trie, long kmer, Integer value) {
					assertNotNull(value);
					assertEquals(value, checkMap.get(kmer));
					checkMap.remove(kmer);
				}
			});
			assertTrue(checkMap.isEmpty());
		} catch (ClassNotFoundException e) {
			throw new RuntimeException(e);
		} catch (IOException e) {
			throw new RuntimeException(e);
		}

	}

	public void testVisit() {
		byte[] cgat = { 'C', 'G', 'A', 'T' };
		Random random = new Random(42);

		KMerStore<Integer> trie = createKMerStore(Integer.class, 35);
		trie.initSize(testSize);
		
		Map<Long, Integer> checkMap = new HashMap<Long, Integer>();
		byte[] read = new byte[trie.getLen()];
		for (int i = 1; i <= testSize; i++) {
			for (int j = 0; j < trie.getLen(); j++) {
				read[j] = cgat[random.nextInt(4)];
			}
			trie.put(read, 0, i, false);

			checkMap.put(CGAT.kmerToLongStraight(read, 0, trie.getLen(), null), i);
		}
		trie.optimize();

		trie.visit(new KMerStoreVisitor<Integer>() {
			@Override
			public void nextValue(KMerStore<Integer> trie, long kmer, Integer value) {
				assertNotNull(value);
				assertEquals(value, checkMap.get(kmer));
				checkMap.remove(kmer);
			}
		});
		assertTrue(checkMap.isEmpty());
	}
}
