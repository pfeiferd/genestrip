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
package org.metagene.genestrip.match;

import java.io.File;
import java.util.Arrays;
import java.util.List;
import java.util.Random;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.atomic.AtomicReference;

import org.junit.Test;
import org.metagene.genestrip.DefaultExecutionContext;
import org.metagene.genestrip.ExecutionContext;
import org.metagene.genestrip.match.FastqKMerMatcher.MatcherReadEntry;
import org.metagene.genestrip.store.Database;
import org.metagene.genestrip.store.KMerSortedArray;
import org.metagene.genestrip.store.KMerStore;
import org.metagene.genestrip.store.KMerStore.ValueConverter;
import org.metagene.genestrip.store.RadixKMerStore;
import org.metagene.genestrip.tax.SmallTaxTree;
import org.metagene.genestrip.tax.SmallTaxTree.SmallTaxIdNode;
import org.metagene.genestrip.tax.TaxTree;
import org.metagene.genestrip.util.CGAT;


import static org.junit.Assert.*;

public class FastqKMerMatcherTest {
	private static final String[] TAXIDS = new String[] { "1", "2", "3" };

	private final Random random = new Random(42);

	// The KMerStore implementations the matcher is tested against (via the common KMerStore interface).
	private enum StoreType {
		SORTED_SMALL, SORTED_LARGE, RADIX
	}

	// Creates an empty KMerStore<String> of the given type, sized to hold exactly the given distinct
	// k-mers, with the taxid values pre-registered.
	private KMerStore<String> newStore(StoreType type, int k, long[] distinctKmers, List<String> initialValues) {
		switch (type) {
		case SORTED_SMALL:
		case SORTED_LARGE: {
			KMerSortedArray<String> store = new KMerSortedArray<>(k, 0.0001, 0.0001, initialValues,
					type == StoreType.SORTED_LARGE, true, distinctKmers.length);
			return store;
		}
		case RADIX: {
			int radixBits = RadixKMerStore.DEFAULT_RADIX_BITS;
			int[] bucketSizes = new int[1 << radixBits];
			for (long kmer : distinctKmers) {
				bucketSizes[RadixKMerStore.radixOf(kmer, radixBits)]++;
			}
			return new RadixKMerStore<>(k, radixBits, bucketSizes, 0.0001, 0.0001, initialValues, true);
		}
		default:
			throw new IllegalArgumentException("Unknown store type: " + type);
		}
	}

	@Test
	public void testMatchRead() {
		// Exercise the matcher against every KMerStore implementation, optimized and not.
		for (StoreType type : StoreType.values()) {
			testMatchReadHelp(type, false);
			testMatchReadHelp(type, true);
		}
	}

	protected void testMatchReadHelp(StoreType type, boolean optimize) {
		int readLength = 500;
		int entries = 2000;

		// Three k-mers in the DB: CC, TT and AG which correspond to the reverse complements GG, AA, CT.
		long ccKmer = CGAT.kMerToLong(new byte[] { 'C', 'C' }, 0, 2, null);
		long ggKmer = CGAT.kMerToLong(new byte[] { 'G', 'G' }, 0, 2, null); // canonicalises to CC
		long ttKmer = CGAT.kMerToLong(new byte[] { 'T', 'T' }, 0, 2, null);
		long agKmer = CGAT.kMerToLong(new byte[] { 'A', 'G' }, 0, 2, null);

		KMerStore<String> store = newStore(type, 2, new long[] { ccKmer, ttKmer, agKmer }, Arrays.asList(TAXIDS));
		store.putLong(ccKmer, TAXIDS[0]);
		assertFalse(store.putLong(ggKmer, TAXIDS[1])); // GG is the reverse complement of CC -> duplicate
		assertTrue(store.putLong(ttKmer, TAXIDS[1]));
		store.putLong(agKmer, TAXIDS[2]);

		if (optimize) {
			store.optimize();
		}

		ExecutionContext bundle = new DefaultExecutionContext(null, 0, 1000);

		FastqKMerMatcher matcher = new MyFastqMatcher(store, readLength * 10, 1, bundle);

		MatcherReadEntry entry = new MatcherReadEntry(2000, true, 4);
		entry.readSize = readLength;

		long[] counters = new long[TAXIDS.length];
		int[] contigs = new int[TAXIDS.length];
		int[] maxContigLen = new int[TAXIDS.length];
		int previousPos = 0;
		int contigLen;
		boolean[] used = new boolean[TAXIDS.length];

		for (int i = 1; i <= entries; i++) {
			// We check correctness of stats for each read separately.
			Arrays.fill(counters, 0);
			Arrays.fill(contigs, 0);
			Arrays.fill(maxContigLen, 0);
			Arrays.fill(used, false);
			matcher.initStats();
			matcher.initUniqueCounter(true);
			entry.bufferPos = 0;
			contigLen = 0;

			int t = -1;
			int lastT = -1;
			byte[] read = entry.read;
			for (int j = 0; j < readLength; j++) {
				read[j] = CGAT.DECODE_TABLE[random.nextInt(4)];
				if (j > 0) {
					lastT = t;
					if ((read[j - 1] == 'C' && read[j] == 'C') || (read[j - 1] == 'G' && read[j] == 'G')) {
						counters[0]++;
						used[0] = true;
						t = 0;
					}
					else if ((read[j - 1] == 'A' && read[j] == 'A') || (read[j - 1] == 'T' && read[j] == 'T')) {
						counters[1]++;
						used[1] = true;
						t = 1;
					}
					else if ((read[j - 1] == 'A' && read[j] == 'G') || (read[j - 1] == 'C' && read[j] == 'T')) {
						counters[2]++;
						used[2] = true;
						t = 2;
					}
					else {
						t = -1;
					}
					if (lastT != t && lastT != -1) {
						contigs[lastT]++;
						if (contigLen > maxContigLen[lastT]) {
							maxContigLen[lastT] = contigLen;
						}
						contigLen = 0;
					}
				}
				if (t != -1) {
					contigLen++;
				}
			}
			if (t != -1) {
				contigs[t]++;
				if (contigLen > maxContigLen[t]) {
					maxContigLen[t] = contigLen;
				}
			}

			/*
			ByteArrayUtil.print(entry.read, System.out);
			System.out.println();
			 */
			matcher.matchRead(entry, 0);
			/*
			System.out.write(entry.buffer, 0, entry.bufferPos);
			System.out.println();
			 */

			matcher.computeUniqueKmerCounts();
			for (int j = 0; j < TAXIDS.length; j++) {
				CountsPerTaxid stats = ((GetStats) matcher).getStats(TAXIDS[j]);
				if (!used[j]) {
					assertNull(stats);
				} else {
					assertEquals(counters[j], stats.getKMers());
					assertEquals(1, stats.getUniqueKMers());
					assertEquals(contigs[j], stats.getContigs());
					assertEquals(maxContigLen[j], stats.getMaxContigLen());
				}
			}
		}
	}

	// Builds the CC/TT/AG store (as in testMatchReadHelp) wrapped in a MyFastqMatcher whose bundle is
	// sized for the given number of consumer threads.
	private MyFastqMatcher newConsistencyMatcher(StoreType type, int initialReadSize, int consumers) {
		long ccKmer = CGAT.kMerToLong(new byte[] { 'C', 'C' }, 0, 2, null);
		long ggKmer = CGAT.kMerToLong(new byte[] { 'G', 'G' }, 0, 2, null); // canonicalises to CC
		long ttKmer = CGAT.kMerToLong(new byte[] { 'T', 'T' }, 0, 2, null);
		long agKmer = CGAT.kMerToLong(new byte[] { 'A', 'G' }, 0, 2, null);

		KMerStore<String> store = newStore(type, 2, new long[] { ccKmer, ttKmer, agKmer }, Arrays.asList(TAXIDS));
		store.putLong(ccKmer, TAXIDS[0]);
		store.putLong(ggKmer, TAXIDS[1]); // duplicate of CC, ignored
		store.putLong(ttKmer, TAXIDS[1]);
		store.putLong(agKmer, TAXIDS[2]);
		store.optimize();

		ExecutionContext bundle = new DefaultExecutionContext(null, consumers, 1000);
		return new MyFastqMatcher(store, initialReadSize, 1, bundle);
	}

	private static void prepareEntry(MatcherReadEntry entry, byte[] read, long readNo) {
		entry.bufferPos = 0;
		entry.usedPaths = 0;
		entry.classNode = null;
		System.arraycopy(read, 0, entry.read, 0, read.length);
		entry.readSize = read.length;
		entry.readNo = readNo;
	}

	// Verifies that matching the same reads across several consumer threads (which share the stats
	// objects) accumulates exactly the same per-tax-id statistics as processing them sequentially.
	// This guards the lock-free/per-contig-batched update path in matchRead() against lost updates.
	@Test
	public void testConcurrentMatchReadConsistency() throws InterruptedException {
		for (StoreType type : StoreType.values()) {
			checkConcurrentConsistency(type);
		}
	}

	private void checkConcurrentConsistency(StoreType type) throws InterruptedException {
		final int readLength = 200;
		final int numReads = 4000;
		final int threads = 8;

		// The same random reads feed both the sequential and the concurrent run.
		Random rnd = new Random(4242);
		byte[][] reads = new byte[numReads][readLength];
		for (int r = 0; r < numReads; r++) {
			for (int j = 0; j < readLength; j++) {
				reads[r][j] = CGAT.DECODE_TABLE[rnd.nextInt(4)];
			}
		}

		// Sequential reference run (single consumer, index 0).
		MyFastqMatcher seq = newConsistencyMatcher(type, readLength * 2, 1);
		seq.initStats();
		MatcherReadEntry seqEntry = new MatcherReadEntry(readLength * 2, true, 4);
		for (int r = 0; r < numReads; r++) {
			prepareEntry(seqEntry, reads[r], r + 1);
			seq.matchRead(seqEntry, 0);
		}

		// Concurrent run: reads partitioned across threads, each with its own consumer index and entry.
		MyFastqMatcher con = newConsistencyMatcher(type, readLength * 2, threads);
		con.initStats();
		CountDownLatch start = new CountDownLatch(1);
		AtomicReference<Throwable> failure = new AtomicReference<>();
		Thread[] ts = new Thread[threads];
		for (int t = 0; t < threads; t++) {
			final int idx = t;
			ts[t] = new Thread(() -> {
				try {
					MatcherReadEntry entry = new MatcherReadEntry(readLength * 2, true, 4);
					start.await();
					for (int r = idx; r < numReads; r += threads) {
						prepareEntry(entry, reads[r], r + 1);
						con.matchRead(entry, idx);
					}
				} catch (Throwable th) {
					failure.compareAndSet(null, th);
				}
			});
			ts[t].start();
		}
		start.countDown();
		for (Thread th : ts) {
			th.join();
		}
		assertNull(failure.get());

		// Aggregate statistics are order-independent sums/maxima, so the concurrent run must match the
		// sequential one exactly.
		for (String tax : TAXIDS) {
			CountsPerTaxid s = seq.getStats(tax);
			CountsPerTaxid c = con.getStats(tax);
			assertNotNull("tax " + tax + " should be hit (" + type + ")", s);
			assertNotNull("tax " + tax + " missing in concurrent run (" + type + ")", c);
			assertEquals("kmers (" + type + ", tax " + tax + ")", s.getKMers(), c.getKMers());
			assertEquals("contigs (" + type + ", tax " + tax + ")", s.getContigs(), c.getContigs());
			assertEquals("reads1KMer (" + type + ", tax " + tax + ")", s.reads1KMer, c.reads1KMer);
			assertEquals("maxContigLen (" + type + ", tax " + tax + ")", s.getMaxContigLen(), c.getMaxContigLen());
		}
	}

	@Test
	public void testReadClassification() {
		for (StoreType type : StoreType.values()) {
			testReadClassificationHelp(type);
		}
	}

	public void testReadClassificationHelp(StoreType type) {
		ClassLoader classLoader = getClass().getClassLoader();
		File treePath = new File(classLoader.getResource("taxtree").getFile());
		TaxTree tree = new TaxTree(treePath, false);
		for (String tax : TAXIDS) {
			tree.getNodeByTaxId(tax).markRequired();
		}
		SmallTaxTree smallTree = tree.toSmallTaxTree();

		// Three k-mers in the DB: CC, CT and CG which correspond to the reverse complements GG, AG, CG.
		long ccKmer = CGAT.kMerToLong(new byte[] { 'C', 'C' }, 0, 2, null);
		long ctKmer = CGAT.kMerToLong(new byte[] { 'C', 'T' }, 0, 2, null);
		long cgKmer = CGAT.kMerToLong(new byte[] { 'C', 'G' }, 0, 2, null);

		KMerStore<String> store = newStore(type, 2, new long[] { ccKmer, ctKmer, cgKmer }, Arrays.asList(TAXIDS));
		store.putLong(ccKmer, TAXIDS[0]);
		assertTrue(store.putLong(ctKmer, TAXIDS[1]));
		store.putLong(cgKmer, TAXIDS[2]);

		Database db = new Database(store, smallTree, null);
		db.initStoreIndices();

		ExecutionContext bundle = new DefaultExecutionContext(null, 0, 1000);
		FastqKMerMatcher matcher = createMatcher2(db.convertKMerStore(), bundle, smallTree, 0);
		matcher.initStats();

		MatcherReadEntry entry = new MatcherReadEntry(10, true, 4);
		entry.readSize = 4;

		fillInRead("CCCC", entry);
		matcher.matchRead(entry, 0);
		assertEquals("1", entry.classNode.getTaxId());
		fillInRead("GAGAGA", entry);
		matcher.matchRead(entry, 0);
		assertNull(entry.classNode);
		fillInRead("CCCG", entry);
		matcher.matchRead(entry, 0);
		assertEquals("3", entry.classNode.getTaxId());
		fillInRead("AGGGG", entry);
		matcher.matchRead(entry, 0);
		assertEquals("2", entry.classNode.getTaxId());
		fillInRead("CCCCCCT", entry);
		matcher.matchRead(entry, 0);
		assertEquals("2", entry.classNode.getTaxId());

		matcher = new MyFastqMatcher2(db.convertKMerStore(), bundle, smallTree, 1);
		matcher.initStats();
		fillInRead("CTCCT", entry);
		matcher.matchRead(entry, 0);
		assertEquals("2", entry.classNode.getTaxId());
		fillInRead("CTCTCCT", entry);
		matcher.matchRead(entry, 0);
		assertNull(entry.classNode);
		fillInRead("TAGGGG", entry);
		matcher.matchRead(entry, 0);
		assertEquals("2", entry.classNode.getTaxId());
		fillInRead("TAGGGGT", entry);
		matcher.matchRead(entry, 0);
		assertNull(entry.classNode);

		// Three k-mers in the DB: CC, CT and CG which correspond to the reverse complements GG, AG, CG.
		matcher = new MyFastqMatcher2(db.convertKMerStore(), bundle, smallTree, 0.5);
		matcher.initStats();
		fillInRead("CCA", entry);
		matcher.matchRead(entry, 0);
		assertEquals("1", entry.classNode.getTaxId());
		fillInRead("CCAA", entry);
		matcher.matchRead(entry, 0);
		assertNull(entry.classNode);

		matcher = new MyFastqMatcher2(db.convertKMerStore(), bundle, smallTree, 0.1);
		matcher.initStats();
		fillInRead("CC", entry);
		matcher.matchRead(entry, 0);
		assertEquals("1", entry.classNode.getTaxId());
		fillInRead("CCA", entry);
		matcher.matchRead(entry, 0);
		assertNull(entry.classNode);
		fillInRead("CCAA", entry);
		matcher.matchRead(entry, 0);
		assertNull(entry.classNode);

		matcher = new MyFastqMatcher2(db.convertKMerStore(), bundle, smallTree, 0.99);
		matcher.initStats();
		fillInRead("TTTT", entry);
		matcher.matchRead(entry, 0);
		assertNull(entry.classNode);
		fillInRead("CTTT", entry);
		matcher.matchRead(entry, 0);
		assertEquals("2", entry.classNode.getTaxId());
	}

	private FastqKMerMatcher createMatcher2(KMerStore<SmallTaxIdNode> kmerStore, ExecutionContext bundle, SmallTaxTree tree,
											double error) {
		return new FastqKMerMatcher (kmerStore, 1024, 100, bundle, false, tree, 4, error, -1, true, 1, null);
	}

	private void fillInRead(String cgat, MatcherReadEntry entry) {
		initEntry(entry);
		for (int i = 0; i < cgat.length(); ++i) {
			entry.read[i] = (byte) cgat.charAt(i);
		}
		entry.read[cgat.length()] = 0;
		entry.readNo++;
		entry.readSize = cgat.length();
	}

	private void initEntry(MatcherReadEntry myEntry) {
		myEntry.bufferPos = 0;
		myEntry.usedPaths = 0;
		myEntry.classNode = null;
		for (int i = 0; i < myEntry.counts.length; i++) {
			myEntry.readTaxIdNode[i] = null;
			myEntry.counts[i] = 0;
		}
	}

	private interface GetStats {
		public CountsPerTaxid getStats(String taxid);
	}

	protected static class MyFastqMatcher extends FastqKMerMatcher implements GetStats {
		private final KMerStore<String> orgKmerStore;

		public MyFastqMatcher(final KMerStore<String> kmerStore, int initialReadSize, int maxQueueSize,
				ExecutionContext bundle) {
			super(kmerStore.convertValues(new ValueConverter<String, SmallTaxIdNode>() {
				@Override
				public SmallTaxIdNode convertValue(String value) {
					SmallTaxIdNode node = new SmallTaxIdNode(value, null, null);
					node.setStoreIndex(kmerStore.getIndexForValue(value));
					return node;
				}
			}), initialReadSize, maxQueueSize, bundle, false, null, 4, 0, 0, true, 1, null);
			orgKmerStore = kmerStore;
			out = System.out;
		}

		@Override
		public CountsPerTaxid getStats(String taxid) {
			return statsIndex[orgKmerStore.getIndexForValue(taxid)];
		}
	}

	protected static class MyFastqMatcher2 extends FastqKMerMatcher {
		public MyFastqMatcher2(KMerStore<SmallTaxIdNode> kmerStore, ExecutionContext bundle, SmallTaxTree tree,
				double error) {
			super(kmerStore, 1024, 100, bundle, false, tree, 4, error, -1, true, 1, null);
		}
	}
}
