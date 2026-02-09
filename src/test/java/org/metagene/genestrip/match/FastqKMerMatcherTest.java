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
import java.util.Random;

import org.junit.Test;
import org.metagene.genestrip.DefaultExecutionContext;
import org.metagene.genestrip.ExecutionContext;
import org.metagene.genestrip.match.FastqKMerMatcher.MatcherReadEntry;
import org.metagene.genestrip.store.Database;
import org.metagene.genestrip.store.KMerSortedArray;
import org.metagene.genestrip.store.KMerSortedArray.ValueConverter;
import org.metagene.genestrip.store.KMerUniqueCounter;
import org.metagene.genestrip.store.KMerUniqueCounterBits;
import org.metagene.genestrip.store.KMerUniqueCounterMap;
import org.metagene.genestrip.tax.SmallTaxTree;
import org.metagene.genestrip.tax.SmallTaxTree.SmallTaxIdNode;
import org.metagene.genestrip.tax.TaxTree;
import org.metagene.genestrip.util.CGAT;

import it.unimi.dsi.fastutil.objects.Object2LongMap;

import static org.junit.Assert.*;

public class FastqKMerMatcherTest {
	private static final String[] TAXIDS = new String[] { "1", "2", "3" };

	private final Random random = new Random(42);

	@Test
	public void testMatchRead() {
		testMatchReadHelp(true,false, false, false);
		testMatchReadHelp(false,false, false, false);
		testMatchReadHelp(true,false, true, false);
		testMatchReadHelp(false,false, true, false);
		testMatchReadHelp(true,false, false, true);
		testMatchReadHelp(false,false, false, true);
		testMatchReadHelp(true,false, true, true);
		testMatchReadHelp(false,false, true, true);
	}

	protected void testMatchReadHelp(boolean inlined, boolean bitMap, boolean large, boolean optimize) {
		int readLength = 500;
		int entries = 2000;

		// Three k-mers in the DB: CC, AA and AG which correspond to the reverse complements GG, TT, CT.
		KMerSortedArray<String> store = new KMerSortedArray<>(2, 0.0001, Arrays.asList(TAXIDS), large, true);
		store.initSize(3);
		byte[] read = new byte[] { 'C', 'C' };
		store.put(read, 0, TAXIDS[0]);
		read[0] = 'G';
		read[1] = 'G';
		assertFalse(store.put(read, 0, TAXIDS[1]));
		read[0] = 'T';
		read[1] = 'T';
		assertTrue(store.put(read, 0, TAXIDS[1]));
		read[0] = 'A';
		read[1] = 'G';
		store.put(read, 0, TAXIDS[2]);

		if (optimize) {
			store.optimize();
		}

		ExecutionContext bundle = new DefaultExecutionContext(null, 0, 1000);

		FastqKMerMatcher matcher = inlined ?
				new MyInlinedFastqMatcher(store, readLength * 10, 1, bundle) :
				new MyFastqMatcher(store, readLength * 10, 1, bundle);
		KMerUniqueCounter uniqueCounter = bitMap ? new KMerUniqueCounterMap() : new KMerUniqueCounterBits(store, true);

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
			matcher.initUniqueCounter(uniqueCounter);
			entry.bufferPos = 0;
			contigLen = 0;

			int t = -1;
			int lastT = -1;
			read = entry.read;
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

			Object2LongMap<String> map = uniqueCounter.getUniqueKmerCounts();
			for (int j = 0; j < TAXIDS.length; j++) {
				CountsPerTaxid stats = ((GetStats) matcher).getStats(TAXIDS[j]);
				if (!used[j]) {
					assertNull(stats);
				} else {
					assertEquals(counters[j], stats.getKMers());
					assertEquals(1, uniqueCounter.getUniqueKmerCount(TAXIDS[j]));
					assertEquals(1, map.getLong(TAXIDS[j]));
					assertEquals(contigs[j], stats.getContigs());
					assertEquals(maxContigLen[j], stats.getMaxContigLen());
				}
			}
		}
	}

	@Test
	public void testReadClassification() {
		testReadClassificationHelp(true, true, true);
		testReadClassificationHelp(false, true, true);
		testReadClassificationHelp(true, false, true);
		testReadClassificationHelp(false, false, true);
		testReadClassificationHelp(true, true, false);
		testReadClassificationHelp(false, true, false);
		testReadClassificationHelp(true, false, false);
		testReadClassificationHelp(false, false, false);
	}

	public void testReadClassificationHelp(boolean inlined, boolean large, boolean optimize) {
		ClassLoader classLoader = getClass().getClassLoader();
		File treePath = new File(classLoader.getResource("taxtree").getFile());
		TaxTree tree = new TaxTree(treePath, false);
		for (String tax : TAXIDS) {
			tree.getNodeByTaxId(tax).markRequired();
		}
		SmallTaxTree smallTree = tree.toSmallTaxTree();

		KMerSortedArray<String> store = new KMerSortedArray<>(2, 0.0001, Arrays.asList(TAXIDS), large, true);
		store.initSize(3);
		byte[] read = new byte[] { 'C', 'C' };
		store.put(read, 0, TAXIDS[0]);
		read[0] = 'C';
		read[1] = 'T';
		assertTrue(store.put(read, 0, TAXIDS[1]));
		read[0] = 'C';
		read[1] = 'G';
		store.put(read, 0, TAXIDS[2]);

		Database db = new Database(store, smallTree);

		ExecutionContext bundle = new DefaultExecutionContext(null, 0, 1000);
		FastqKMerMatcher matcher = createMatcher2(inlined, db.convertKMerStore(), bundle, smallTree, 0);
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

	private FastqKMerMatcher createMatcher2(boolean inlined, KMerSortedArray<SmallTaxIdNode> kmerStore, ExecutionContext bundle, SmallTaxTree tree,
											double error) {
		if (inlined) {
			return new InlinedFastqKMerMatcher (kmerStore, 1024, 100, bundle, false, 10, tree, 4, error, -1, true);
		}
		else {
			return new FastqKMerMatcher (kmerStore, 1024, 100, bundle, false, 10, tree, 4, error, -1, true);
		}
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
		private final KMerSortedArray<String> orgKmerStore;

		public MyFastqMatcher(KMerSortedArray<String> kmerStore, int initialReadSize, int maxQueueSize,
				ExecutionContext bundle) {
			super(new KMerSortedArray<>(kmerStore, new ValueConverter<String, SmallTaxIdNode>() {
				@Override
				public SmallTaxIdNode convertValue(String value) {
					SmallTaxIdNode node = new SmallTaxIdNode(value, null, null);
					node.setStoreIndex(kmerStore.getIndexForValue(value));
					return node;
				}
			}), initialReadSize, maxQueueSize, bundle, false, 10, null, 4, 0, 0, true);
			orgKmerStore = kmerStore;
			out = System.out;
		}

		@Override
		public CountsPerTaxid getStats(String taxid) {
			return statsIndex[orgKmerStore.getIndexForValue(taxid)];
		}
	}

	protected static class MyInlinedFastqMatcher extends InlinedFastqKMerMatcher implements GetStats {
		private final KMerSortedArray<String> orgKmerStore;

		public MyInlinedFastqMatcher(KMerSortedArray<String> kmerStore, int initialReadSize, int maxQueueSize,
							  ExecutionContext bundle) {
			super(new KMerSortedArray<>(kmerStore, new ValueConverter<String, SmallTaxIdNode>() {
				@Override
				public SmallTaxIdNode convertValue(String value) {
					SmallTaxIdNode node = new SmallTaxIdNode(value, null, null);
					node.setStoreIndex(kmerStore.getIndexForValue(value));
					return node;
				}
			}), initialReadSize, maxQueueSize, bundle, false, 10, null, 4, 0, 0, true);
			orgKmerStore = kmerStore;
			out = System.out;
		}

		@Override
		public CountsPerTaxid getStats(String taxid) {
			return statsIndex[orgKmerStore.getIndexForValue(taxid)];
		}
	}

	protected static class MyFastqMatcher2 extends FastqKMerMatcher {
		public MyFastqMatcher2(KMerSortedArray<SmallTaxIdNode> kmerStore, ExecutionContext bundle, SmallTaxTree tree,
				double error) {
			super(kmerStore, 1024, 100, bundle, false, 10, tree, 4, error, -1, true);
		}
	}

	protected static class MyInlinedFastqMatcher2 extends InlinedFastqKMerMatcher {
		public MyInlinedFastqMatcher2(KMerSortedArray<SmallTaxIdNode> kmerStore, ExecutionContext bundle, SmallTaxTree tree,
							   double error) {
			super(kmerStore, 1024, 100, bundle, false, 10, tree, 4, error, -1, true);
		}
	}
}
