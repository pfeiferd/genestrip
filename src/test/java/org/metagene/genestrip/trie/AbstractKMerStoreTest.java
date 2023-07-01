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
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
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
import org.metagene.genestrip.util.CGATRingBuffer;
import org.metagene.genestrip.util.StreamProvider;

import junit.framework.TestCase;

public abstract class AbstractKMerStoreTest extends TestCase implements KMerStoreFactory {
	protected final int k = 31;
	protected final int testSize = 5 * 1000 * 1000;
	protected final int negativeTestSize = testSize;
	protected final Random random = new Random(42);

	public void testChain() throws IOException {
		Main main = new Main();
		main.parseAndRun(new String[] { "bart_h", "clear", "genall" });

		GSProject project = main.getProject();

		File bartHReads = ((FileGoal<GSProject>) main.getMaker().getGoal("kmerfastq")).getFile();

		String outGoal = project.getConfig().isUseKraken1() ? "sort" : "kmerkrakenout";
		File fromKraken = ((FileGoal<GSProject>) main.getMaker().getGoal(outGoal)).getFile();

		final KMerStore<String> store = createKMerStore(String.class, project.getkMserSize());

		MurmurCGATBloomFilter bloomFilter = new MurmurCGATBloomFilter(store.getLen(), 0.00001);
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
						store.put(read, 0, kmerTaxid, false);
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
		if (store instanceof KMerTrie) {
			checkTrie(fromKraken, bartHReads, krakenFilter, store, bloomFilter, nodes);
		}

		store.optimize();

		// Test compressed;
		checkTrie(fromKraken, bartHReads, krakenFilter, store, bloomFilter, nodes);
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
		for (int i = 1; i <= negativeTestSize; i++) {
			for (int j = 0; j < read.length; j++) {
				read[j] = CGAT.DECODE_TABLE[random.nextInt(4)];
			}
			if (!bloomFilter.containsStraight(read, 0, null)) {
				assertNull(trie.get(read, 0, false));
			}
		}
	}

	public void testSaveLoad() {
		KMerStore<Integer> store = createKMerStore(Integer.class, k);
		store.initSize(testSize);
		
		Map<List<Byte>, Integer> controlMap = new HashMap<List<Byte>, Integer>();
		Map<Long, Integer> kmerMap = new HashMap<Long, Integer>();
		fillStore(testSize, store, controlMap, kmerMap);
		store.optimize();
		
		try {
			// Check serialization:
			File tmpdir = Files.createTempDirectory(Paths.get("target"), "serialization-test").toFile();
			tmpdir.deleteOnExit();
			File saved = new File(tmpdir, "saved_trie.ser");
			saved.deleteOnExit();
			if (saved.exists()) {
				saved.delete();
			}
			KMerStore.save(store, saved);

			KMerStore<Integer> loadedStore = KMerStore.load(saved);
			checkStoreContent(testSize, loadedStore, controlMap);
			checkVisitation(loadedStore, kmerMap);			
		} catch (ClassNotFoundException e) {
			throw new RuntimeException(e);
		} catch (IOException e) {
			throw new RuntimeException(e);
		}

	}

	protected void fillStore(int nValues, KMerStore<Integer> store, Map<List<Byte>, Integer> controlMap,
			Map<Long, Integer> controlMap2) {
		byte[] read = new byte[store.getLen()];
		List<Byte> readAsList = null;
		for (int i = 1; i <= nValues; i++) {
			if (controlMap != null) {
				readAsList = new ArrayList<Byte>();
			}
			for (int j = 0; j < store.getLen(); j++) {
				read[j] = CGAT.DECODE_TABLE[random.nextInt(4)];
				if (controlMap != null) {
					readAsList.add(read[j]);
				}
			}
			store.put(read, 0, i, false);
			if (controlMap != null) {
				controlMap.put(readAsList, i);
			}
			if (controlMap2 != null) {
				controlMap2.put(CGAT.kmerToLongStraight(read, 0, k, null), i);
			}
		}
	}

	protected void checkStoreContent(int negTestValues, KMerStore<Integer> store, Map<List<Byte>, Integer> controlMap) {
		byte[] read = new byte[store.getLen()];
		CGATRingBuffer ringBuffer = new CGATRingBuffer(read.length);

		for (List<Byte> key : controlMap.keySet()) {
			Integer v = controlMap.get(key);
			for (int j = 0; j < read.length; j++) {
				read[j] = CGAT.DECODE_TABLE[random.nextInt(4)];
				ringBuffer.put(read[j]);
			}

			assertEquals(v, store.get(read, 0, false));
			assertEquals(v, store.get(ringBuffer, false));

			CGAT.reverse(read);
			for (int j = 0; j < read.length; j++) {
				ringBuffer.put(read[j]);
			}
			assertEquals(v, store.get(read, 0, true));
			assertEquals(v, store.get(ringBuffer, true));
		}

		List<Byte> readAsList = new ArrayList<Byte>();
		for (int i = 1; i <= negativeTestSize; i++) {
			readAsList.clear();
			for (int j = 0; j < read.length; j++) {
				read[j] = CGAT.DECODE_TABLE[random.nextInt(4)];
				readAsList.add(read[j]);
				ringBuffer.put(read[j]);
			}
			if (!controlMap.containsKey(readAsList)) {
				assertNull(store.get(read, 0, false));
				assertNull(store.get(ringBuffer, false));
			}

			CGAT.reverse(read);
			readAsList.clear();
			for (int j = 0; j < read.length; j++) {
				readAsList.add(read[j]);
				ringBuffer.put(read[j]);
			}
			if (!controlMap.containsKey(readAsList)) {
				assertNull(store.get(read, 0, true));
				assertNull(store.get(ringBuffer, true));
			}
		}
	}
	
	public void testPutGet() {
		KMerStore<Integer> store = createKMerStore(Integer.class, k);
		store.initSize(testSize);
		
		Map<List<Byte>, Integer> controlMap = new HashMap<List<Byte>, Integer>();
		fillStore(testSize, store, controlMap, null);
		store.optimize();
		
		checkStoreContent(testSize, store, controlMap);
	}

	public void testVisit() {
		KMerStore<Integer> store = createKMerStore(Integer.class, k);
		store.initSize(testSize);
		
		Map<Long, Integer> kmerMap = new HashMap<Long, Integer>();
		fillStore(testSize, store, null, kmerMap);
		store.optimize();
		
		checkVisitation(store, kmerMap);
	}
	
	protected void checkVisitation(KMerStore<Integer> store, Map<Long, Integer> kmerMap) {
		store.visit(new KMerStoreVisitor<Integer>() {
			@Override
			public void nextValue(KMerStore<Integer> trie, long kmer, Integer value) {
				assertNotNull(value);
				assertEquals(value, kmerMap.get(kmer));
				kmerMap.remove(kmer);
			}
		});
		assertTrue(kmerMap.isEmpty());		
	}
}
