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

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;

import java.io.File;
import java.util.Arrays;
import java.util.Random;

import org.junit.Test;
import org.metagene.genestrip.DefaultExecutionContext;
import org.metagene.genestrip.ExecutionContext;
import org.metagene.genestrip.match.FastqKMerMatcher.MyReadEntry;
import org.metagene.genestrip.store.Database;
import org.metagene.genestrip.store.KMerSortedArray;
import org.metagene.genestrip.store.KMerSortedArray.ValueConverter;
import org.metagene.genestrip.store.KMerUniqueCounter;
import org.metagene.genestrip.store.KMerUniqueCounterBits;
import org.metagene.genestrip.store.KMerUniqueCounterMap;
import org.metagene.genestrip.tax.SmallTaxTree;
import org.metagene.genestrip.tax.SmallTaxTree.SmallTaxIdNode;
import org.metagene.genestrip.tax.TaxTree;
import org.metagene.genestrip.util.ByteArrayUtil;
import org.metagene.genestrip.util.CGAT;

import it.unimi.dsi.fastutil.objects.Object2LongMap;

public class FastqKMerMatcherTest {
	private static final String[] TAXIDS = new String[] { "1", "2", "3" };

	private final Random random = new Random(42);

	@Test
	public void testMatchRead() {
		//testMatchReadHelp(true);
		testMatchReadHelp(false);
	}

	protected void testMatchReadHelp(boolean bitMap) {
		int readLength = 5;
		int entries = 2000;

		KMerSortedArray<String> store = new KMerSortedArray<>(1, 0.0001, Arrays.asList(TAXIDS), false);
		store.initSize(3);
		byte[] read = new byte[] { 'C' };
		store.put(read, 0, TAXIDS[0], false);
		read[0] = 'G';
		store.put(read, 0, TAXIDS[1], false);
		read[0] = 'A';
		store.put(read, 0, TAXIDS[2], false);

		ExecutionContext bundle = new DefaultExecutionContext(0, 1000);

		MyFastqMatcher matcher = new MyFastqMatcher(store, readLength * 10, 1, bundle);
		KMerUniqueCounter uniqueCounter = bitMap ? new KMerUniqueCounterMap() : new KMerUniqueCounterBits(store, true);

		MyReadEntry entry = new MyReadEntry(2000, true, 4);
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
			matcher.initRoot();
			matcher.initUniqueCounter(uniqueCounter);
			entry.bufferPos = 0;
			contigLen = 0;

			for (int j = 0; j < readLength; j++) {
				int pos = random.nextInt(4);
				entry.read[j] = CGAT.DECODE_TABLE[pos];
				if (pos != 3) {
					counters[pos]++;
					used[pos] = true;
				}
				if (j > 0) {
					if (pos != previousPos) {
						if (previousPos != 3) {
							contigs[previousPos]++;
							if (contigLen > maxContigLen[previousPos]) {
								maxContigLen[previousPos] = contigLen;
							}
						}
						contigLen = 0;
					}
				}
				contigLen++;
				previousPos = pos;
			}
			if (previousPos != 3) {
				contigs[previousPos]++;
				if (contigLen > maxContigLen[previousPos]) {
					maxContigLen[previousPos] = contigLen;
				}
			}
			matcher.matchRead(entry, i, false);

			ByteArrayUtil.print(entry.read, System.out);
			System.out.println();
			System.out.write(entry.buffer, 0, entry.bufferPos);
			System.out.println();

			Object2LongMap<String> map = uniqueCounter.getUniqueKmerCounts();
			for (int j = 0; j < TAXIDS.length; j++) {
				CountsPerTaxid stats = matcher.getStats(TAXIDS[j]);
				if (!used[j]) {
					assertNull(stats);
				} else {
					assertEquals(counters[j], stats.getKMers());
					assertEquals(used[j] ? 1 : 0, uniqueCounter.getUniqueKmerCount(TAXIDS[j]));
					assertEquals(used[j] ? 1 : 0, map.getLong(TAXIDS[j]));
					assertEquals(contigs[j], stats.getContigs());
					assertEquals(maxContigLen[j], stats.getMaxContigLen());
				}
			}
		}
	}

	@Test
	public void testReadClassification() {
		ClassLoader classLoader = getClass().getClassLoader();
		File treePath = new File(classLoader.getResource("taxtree").getFile());
		TaxTree tree = new TaxTree(treePath);
		for (String tax : TAXIDS) {
			tree.getNodeByTaxId(tax).markRequired();
		}
		SmallTaxTree smallTree = tree.toSmallTaxTree();

		KMerSortedArray<String> store = new KMerSortedArray<>(1, 0.0001, Arrays.asList(TAXIDS), false);
		store.initSize(3);
		byte[] read = new byte[] { 'C' };
		store.put(read, 0, TAXIDS[0], false);
		read[0] = 'G';
		store.put(read, 0, TAXIDS[1], false);
		read[0] = 'A';
		store.put(read, 0, TAXIDS[2], false);
		store.optimize();

		Database db = new Database(store, smallTree);

		ExecutionContext bundle = new DefaultExecutionContext(0, 1000);
		MyFastqMatcher2 matcher = new MyFastqMatcher2(db.convertKMerStore(), bundle, smallTree, 0);
		matcher.initRoot();

		MyReadEntry entry = new MyReadEntry(10, true, 4);
		entry.readSize = 4;

		fillInRead("TTTT", entry);
		matcher.matchRead(entry, 0, false);
		assertNull(entry.classNode);
		fillInRead("CCCT", entry);
		matcher.matchRead(entry, 0, false);
		assertNull(entry.classNode);
		fillInRead("CCCC", entry);
		matcher.matchRead(entry, 0, false);
		assertEquals("1", entry.classNode.getTaxId());
		fillInRead("CCCG", entry);
		matcher.matchRead(entry, 0, false);
		assertEquals("2", entry.classNode.getTaxId());
		fillInRead("CGGG", entry);
		matcher.matchRead(entry, 0, false);
		assertEquals("2", entry.classNode.getTaxId());
		fillInRead("CAAA", entry);
		matcher.matchRead(entry, 0, false);
		assertEquals("3", entry.classNode.getTaxId());
		fillInRead("GCAG", entry);
		matcher.matchRead(entry, 0, false);
		assertEquals("2", entry.classNode.getTaxId());
		fillInRead("ACAG", entry);
		matcher.matchRead(entry, 0, false);
		assertEquals("3", entry.classNode.getTaxId());
		fillInRead("GCAA", entry);
		matcher.matchRead(entry, 0, false);
		assertEquals("3", entry.classNode.getTaxId());
		fillInRead("CCAG", entry);
		matcher.matchRead(entry, 0, false);
		assertEquals("1", entry.classNode.getTaxId());

		matcher = new MyFastqMatcher2(db.convertKMerStore(), bundle, smallTree, 1);
		matcher.initRoot();
		fillInRead("CCCT", entry);
		matcher.matchRead(entry, 0, false);
		assertEquals("1", entry.classNode.getTaxId());
		fillInRead("TCCT", entry);
		matcher.matchRead(entry, 0, false);
		assertNull(entry.classNode);

		matcher = new MyFastqMatcher2(db.convertKMerStore(), bundle, smallTree, 0.5);
		matcher.initRoot();
		fillInRead("CCCT", entry);
		matcher.matchRead(entry, 0, false);
		assertEquals("1", entry.classNode.getTaxId());
		fillInRead("TCCT", entry);
		matcher.matchRead(entry, 0, false);
		assertEquals("1", entry.classNode.getTaxId());
		fillInRead("TTCT", entry);
		matcher.matchRead(entry, 0, false);
		assertNull(entry.classNode);

		matcher = new MyFastqMatcher2(db.convertKMerStore(), bundle, smallTree, 0.1);
		matcher.initRoot();
		fillInRead("CCCT", entry);
		matcher.matchRead(entry, 0, false);
		assertNull(entry.classNode);
		fillInRead("CCCC", entry);
		matcher.matchRead(entry, 0, false);
		assertEquals("1", entry.classNode.getTaxId());

		matcher = new MyFastqMatcher2(db.convertKMerStore(), bundle, smallTree, 0.99);
		matcher.initRoot();
		fillInRead("TTTT", entry);
		matcher.matchRead(entry, 0, false);
		assertNull(entry.classNode);
		fillInRead("CTTT", entry);
		matcher.matchRead(entry, 0, false);
		assertEquals("1", entry.classNode.getTaxId());
	}

	private void fillInRead(String cgat, MyReadEntry entry) {
		initEntry(entry);
		for (int i = 0; i < cgat.length(); ++i) {
			entry.read[i] = (byte) cgat.charAt(i);
		}
		entry.read[cgat.length()] = 0;
		entry.readNo++;
	}

	private void initEntry(MyReadEntry myEntry) {
		myEntry.bufferPos = 0;
		myEntry.usedPaths = 0;
		myEntry.classNode = null;
		for (int i = 0; i < myEntry.counts.length; i++) {
			myEntry.readTaxIdNode[i] = null;
			myEntry.counts[i] = 0;
		}
	}

	protected static class MyFastqMatcher extends FastqKMerMatcher {
		private final KMerSortedArray<String> orgKmerStore;

		public MyFastqMatcher(KMerSortedArray<String> kmerStore, int initialReadSize, int maxQueueSize,
				ExecutionContext bundle) {
			super(new KMerSortedArray<>(kmerStore, new ValueConverter<String, SmallTaxIdNode>() {
				@Override
				public SmallTaxIdNode convertValue(String value) {
					SmallTaxIdNode node = new SmallTaxIdNode(value);
					node.setStoreIndex(kmerStore.getIndexForValue(value));
					return node;
				}
			}), initialReadSize, maxQueueSize, bundle, false, 10, null, 4, 0);
			orgKmerStore = kmerStore;
			out = System.out;
		}

		@Override
		public void initRoot() {
			super.initRoot();
		}

		@Override
		public void initUniqueCounter(KMerUniqueCounter uniqueCounter) {
			super.initUniqueCounter(uniqueCounter);
		}

		public CountsPerTaxid getStats(String taxid) {
			return statsIndex[orgKmerStore.getIndexForValue(taxid)];
		}

		/*
		@Override
		protected boolean matchRead(MyReadEntry entry, int index, boolean reverse) {
			return super.matchRead(entry, index, reverse);
		}
		*/
	}

	protected static class MyFastqMatcher2 extends FastqKMerMatcher {
		public MyFastqMatcher2(KMerSortedArray<SmallTaxIdNode> kmerStore, ExecutionContext bundle, SmallTaxTree tree,
				double error) {
			super(kmerStore, 1024, 100, bundle, false, 10, tree, 4, error);
		}

		@Override
		public void initRoot() {
			super.initRoot();
		}

		/*
		@Override
		protected boolean matchRead(MyReadEntry entry, int index, boolean reverse) {
			return super.matchRead(entry, index, reverse);
		}
		*/
	}
}
