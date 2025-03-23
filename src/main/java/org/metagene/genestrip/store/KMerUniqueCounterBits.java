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

import java.util.HashMap;
import java.util.Map;

import org.metagene.genestrip.store.KMerSortedArray.KMerSortedArrayVisitor;
import org.metagene.genestrip.util.LargeBitVector;
import org.metagene.genestrip.util.LargeShortVector;

import it.unimi.dsi.fastutil.objects.Object2LongLinkedOpenHashMap;
import it.unimi.dsi.fastutil.objects.Object2LongMap;

public class KMerUniqueCounterBits implements KMerUniqueCounter {
	private final KMerSortedArray<String> store;
	private final LargeBitVector bitVector;
	private final LargeShortVector countsVector;

	public KMerUniqueCounterBits(KMerSortedArray<String> store, boolean withCounts) {
		this.store = store;
		bitVector = new LargeBitVector(store.getEntries());
		countsVector = withCounts ? new LargeShortVector(store.getEntries()) : null;
	}
	
	public boolean isWithCounts() {
		return countsVector != null;
	}

	@Override
	public void clear() {
		bitVector.clear();
		if (countsVector != null) {
			countsVector.clear();
		}
	}

	@Override
	public final synchronized void put(final long kmer, final String taxid, final long index) {
		bitVector.set(index);
		if (countsVector != null) {
			countsVector.inc(index);
		}
	}

	public final synchronized void putInlined(final long kmer, final String taxid, final long index) {
		if (bitVector.largeBits != null) {
			long arrayIndex = ((index >>> 6) % bitVector.size);
			bitVector.largeBits[(int) (arrayIndex >>> 27)][(int) (arrayIndex & 134217727)] = bitVector.largeBits[(int) (arrayIndex >>> 27)][(int) (arrayIndex & 134217727)] | (1L << (index & 0b111111));
		} else {
			int arrayIndex = (int) ((index >>> 6) % bitVector.size);
			bitVector.bits[arrayIndex] = bitVector.bits[arrayIndex] | (1L << (index & 0b111111));
		}
		if (countsVector != null) {
			if (countsVector.largeShorts != null) {
				countsVector.largeShorts[(int) (index >>> 27)][(int) (index & 134217727)]++;
			} else {
				++countsVector.shorts[(int) index];
			}
		}
	}

	@Override
	public Object2LongMap<String> getUniqueKmerCounts() {
		Object2LongMap<String> res = new Object2LongLinkedOpenHashMap<String>();
		long[] valueCounter = new long[store.getNValues()];
		store.visit(new KMerSortedArrayVisitor<String>() {
			@Override
			public void nextValue(KMerSortedArray<String> trie, long kmer, short index, long i) {
				if (bitVector.get(i)) {
					valueCounter[index]++;
				}
			}
		});
		for (short i = 0; i < valueCounter.length; i++) {
			String taxid = store.getValueForIndex((short)i);
			res.put(taxid, valueCounter[i]);
		}

		return res;
	}

	public Map<String, short[]> getMaxCountsCounts(int counts) {
		Map<String, short[]> res = new HashMap<String, short[]>();
		short[] totalMaxCounts = new short[counts];
		res.put(null, totalMaxCounts);
		
		
		store.visit(new KMerSortedArrayVisitor<String>() {
			@Override
			public void nextValue(KMerSortedArray<String> trie, long kmer, short index, long i) {
				if (bitVector.get(i)) {
					String taxid = store.getValueForIndex(index);
					if (taxid != null) {						
						short[] target = res.get(taxid); 
						if (target == null) {
							target = new short[counts];
							res.put(taxid, target);
						}								
						short count = countsVector.get(i);
						updateMaxCounts(count, target);
						updateMaxCounts(count, totalMaxCounts);
					}
				}
			}
		});

		return res;
	}
	
	private void updateMaxCounts(short count, short[] target) {
		for (int j = 0; j < target.length; j++) {
			if (count > target[j]) {
				for (int k = target.length - 1; k > j; k--) {
					target[k] = target[k - 1];
				}
				target[j] = count;
				return;
			}
		}		
	}

	@Override
	public int getUniqueKmerCount(String taxid) {
		final short sindex = store.getIndexForValue(taxid);
		if (sindex < 0) {
			return 0;
		}
		int[] count = new int[1];
		store.visit(new KMerSortedArrayVisitor<String>() {
			@Override
			public void nextValue(KMerSortedArray<String> trie, long kmer, short index, long i) {
				if (index == sindex && bitVector.get(i)) {
					count[0]++;
				}
			}
		});

		return count[0];
	}

	public void getMaxCounts(String taxid, short[] target) {
		final short sindex = store.getIndexForValue(taxid);
		if (sindex < 0) {
			return;
		}
		store.visit(new KMerSortedArrayVisitor<String>() {
			@Override
			public void nextValue(KMerSortedArray<String> trie, long kmer, short index, long i) {
				if (index == sindex && bitVector.get(i)) {
					updateMaxCounts(countsVector.get(i), target);
				}
			}
		});
	}
}