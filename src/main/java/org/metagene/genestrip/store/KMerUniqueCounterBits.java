package org.metagene.genestrip.store;

import org.metagene.genestrip.store.KMerSortedArray.KMerSortedArrayVisitor;
import org.metagene.genestrip.util.LargeBitVector;

import it.unimi.dsi.fastutil.objects.Object2LongLinkedOpenHashMap;
import it.unimi.dsi.fastutil.objects.Object2LongMap;

public class KMerUniqueCounterBits implements KMerUniqueCounter {
	private final KMerSortedArray<String> store;
	private final LargeBitVector bitVector;

	public KMerUniqueCounterBits(KMerSortedArray<String> store) {
		this.store = store;
		bitVector = new LargeBitVector(store.getEntries());
	}

	@Override
	public void clear() {
		bitVector.clear();
	}

	@Override
	public synchronized void put(long kmer, String taxid, long index) {
		bitVector.set(index);
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
}