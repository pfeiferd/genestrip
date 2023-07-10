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
package org.metagene.genestrip.store;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

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

import org.junit.Ignore;
import org.junit.Test;
import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.Main;
import org.metagene.genestrip.bloom.MurmurCGATBloomFilter;
import org.metagene.genestrip.goals.KrakenResCountGoal;
import org.metagene.genestrip.kraken.KrakenResultFastqMergeListener;
import org.metagene.genestrip.kraken.KrakenResultFastqMerger;
import org.metagene.genestrip.make.FileGoal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.store.KMerStore.KMerStoreVisitor;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;
import org.metagene.genestrip.util.CGAT;
import org.metagene.genestrip.util.CGATRingBuffer;
import org.metagene.genestrip.util.StreamProvider;

public abstract class AbstractKMerStoreTest implements KMerStoreFactory {
	protected final Random random = new Random(42);

	protected int k = 31;
	protected int testSize = 1 * 1000 * 1000;
	protected int negativeTestSize = testSize;

	@Ignore
	@Test
	public void testChain() throws IOException {
		Main main = new Main();
		main.parseAndRun(new String[] { "bart_h", "clear", "genall" });

		GSProject project = main.getProject();

		File bartHReads = ((FileGoal<GSProject>) main.getMaker().getGoal("kmerfastq")).getFile();

		String outGoal = project.getConfig().isUseKraken1() ? "sort" : "kmerkrakenout";
		File fromKraken = ((FileGoal<GSProject>) main.getMaker().getGoal(outGoal)).getFile();

		final KMerStore<String> store = createKMerStore(String.class, project.getKMserSize());

		MurmurCGATBloomFilter bloomFilter = new MurmurCGATBloomFilter(store.getK(), 0.00001);
		bloomFilter.clear();
		bloomFilter.ensureExpectedSize(5 * 1000 * 1000, false);

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
		byte[] read = new byte[trie.getK()];
		for (int i = 1; i <= negativeTestSize; i++) {
			for (int j = 0; j < read.length; j++) {
				read[j] = CGAT.DECODE_TABLE[random.nextInt(4)];
			}
			if (!bloomFilter.contains(read, 0, null, false)) {
				assertNull(trie.get(read, 0, false));
			}
		}
	}

	@Test
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
			if (store instanceof KMerSortedArray) {
				checkVisitation(loadedStore, kmerMap);
			}
		} catch (ClassNotFoundException e) {
			throw new RuntimeException(e);
		} catch (IOException e) {
			throw new RuntimeException(e);
		}

	}

	protected void fillStore(int nEntries, KMerStore<Integer> store, Map<List<Byte>, Integer> controlMap,
			Map<Long, Integer> controlMap2) {
		byte[] read = new byte[store.getK()];
		List<Byte> readAsList = null;
		for (int i = 0; i < nEntries; i++) {
			if (controlMap != null) {
				readAsList = new ArrayList<Byte>();
			}
			for (int j = 0; j < store.getK(); j++) {
				read[j] = CGAT.DECODE_TABLE[random.nextInt(4)];
				if (controlMap != null) {
					readAsList.add(read[j]);
				}
			}
			int v = i % (KMerSortedArray.MAX_VALUES - 1);
			store.put(read, 0, v, false);
			if (controlMap != null) {
				controlMap.put(readAsList, v);
			}
			if (controlMap2 != null) {
				controlMap2.put(CGAT.kMerToLongStraight(read, 0, k, null), v);
			}
		}
	}

	protected void checkStoreContent(int negTestValues, KMerStore<Integer> store, Map<List<Byte>, Integer> controlMap) {
		byte[] read = new byte[store.getK()];
		CGATRingBuffer ringBuffer = new CGATRingBuffer(read.length);

		for (List<Byte> key : controlMap.keySet()) {
			Integer v = controlMap.get(key);
			for (int j = 0; j < read.length; j++) {
				read[j] = key.get(j);
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

	@Test
	public void testPutGet() {
		KMerStore<Integer> store = createKMerStore(Integer.class, k);
		store.initSize(testSize);

		Map<List<Byte>, Integer> controlMap = new HashMap<List<Byte>, Integer>();
		fillStore(testSize, store, controlMap, null);
		store.optimize();

		checkStoreContent(testSize, store, controlMap);
	}

	@Test
	public void testVisit() {
		KMerStore<Integer> store = createKMerStore(Integer.class, k);
		store.initSize(testSize);

		Map<Long, Integer> kmerMap = new HashMap<Long, Integer>();
		fillStore(testSize, store, null, kmerMap);
		if (store instanceof KMerSortedArray) {
			store.optimize();
		}

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
