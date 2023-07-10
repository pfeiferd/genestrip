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

import java.util.Arrays;
import java.util.Random;

import org.junit.Test;
import org.metagene.genestrip.match.FastqKMerMatcher.MyReadEntry;
import org.metagene.genestrip.match.ResultReporter.StatsPerTaxid;
import org.metagene.genestrip.store.KMerSortedArray;
import org.metagene.genestrip.store.KMerUniqueCounterBits;
import org.metagene.genestrip.util.ByteArrayUtil;
import org.metagene.genestrip.util.CGAT;

public class FastqKMerMatcherTest {
	private static final String[] TAXIDS = new String[] { "1", "2", "3" };
	
	private Random random = new Random(42);

	@Test
	public void testMatchRead() {
		int readLength = 150;
		int entries = 200;
		
		KMerSortedArray<String> store = new KMerSortedArray<String>(1, 0.0001, Arrays.asList(TAXIDS), false);
		byte[] read = new byte[] { 'C' };
		store.put(read, 0, TAXIDS[0], false);
		read[0] = 'G';
		store.put(read, 0, TAXIDS[0], false);
		read[0] = 'A';
		store.put(read, 0, TAXIDS[0], false);

		MyFastqMatcher matcher = new MyFastqMatcher(store, readLength, 1, 0);
		KMerUniqueCounterBits uniqueCounterMap = new KMerUniqueCounterBits(store);

		MyReadEntry entry = new MyReadEntry(2000, 1);
		entry.readSize = readLength;

		long[] counters = new long[TAXIDS.length];
		int[] contigs = new int[TAXIDS.length];
		int[] maxContigLen = new int[TAXIDS.length];
		int previousPos = 0;
		int contigLen;

		for (int i = 1; i <= entries; i++) {
			// We check correctness of stats for each read separately.
			Arrays.fill(counters, 0);
			Arrays.fill(contigs, 0);
			Arrays.fill(maxContigLen, 0);
			matcher.initRoot();
			uniqueCounterMap.clear();
			entry.bufferPos = 0;
			contigLen = 0;

			for (int j = 0; j < readLength; j++) {
				int pos = random.nextInt(4);
				entry.read[j] = CGAT.DECODE_TABLE[pos];
				System.out.print((char) entry.read[j]);
				if (pos != 3) {
					counters[pos]++;
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
			ByteArrayUtil.print(entry.read, System.out);
			matcher.classifyRead(entry, false);

			ByteArrayUtil.print(entry.read, System.out);
			System.out.println();
			System.out.write(entry.buffer, 0, entry.bufferPos);
			System.out.println();

			for (int j = 0; j < TAXIDS.length; j++) {
				StatsPerTaxid stats = matcher.getStats(TAXIDS[j]);
				assertEquals(stats.getKMers(), counters[j]);
				assertEquals(counters[j], uniqueCounterMap.getUniqueKmerCount(TAXIDS[j]));
				assertEquals(contigs[j], stats.getContigs());
				assertEquals(maxContigLen[j], stats.getMaxContigLen());
			}
		}
	}

	protected static class MyFastqMatcher extends FastqKMerMatcher {
		public MyFastqMatcher(KMerSortedArray<String> kmerStore, int maxReadSize, int maxQueueSize, int consumerNumber) {
			super(kmerStore, maxReadSize, maxQueueSize, consumerNumber);
			out = System.out;
		}

		@Override
		public void initRoot() {
			super.initRoot();
		}

		public StatsPerTaxid getStats(String taxid) {
			return root.get(taxid);
		}

		@Override
		public boolean classifyRead(MyReadEntry entry, boolean reverse) {
			return super.classifyRead(entry, reverse);
		}
	}
}
