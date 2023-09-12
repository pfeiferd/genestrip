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
	public synchronized void put(long kmer, String taxid, long index) {
		bitVector.set(index);
		if (countsVector != null) {
			countsVector.inc(index);
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
			String taxid = store.getValueForIndex(i);
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
				for (int k = j + 1; k < target.length; k++) {
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