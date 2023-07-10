package org.metagene.genestrip.store;

import org.metagene.genestrip.store.KMerSortedArray.KMerSortedArrayVisitor;
import org.metagene.genestrip.util.LargeBitVector;

import it.unimi.dsi.fastutil.objects.Object2IntLinkedOpenHashMap;
import it.unimi.dsi.fastutil.objects.Object2IntMap;

public class KMerUniqueCounterBits implements KMerUniqueCounter {
	private final KMerSortedArray<String> store;
	private final LargeBitVector bitVector;

	public KMerUniqueCounterBits(KMerSortedArray<String> store) {
		this.store = store;
		bitVector = new LargeBitVector(store.getSize());
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
	public Object2IntMap<String> getUniqueKmerCounts() {
		Object2IntMap<String> res = new Object2IntLinkedOpenHashMap<String>();
		int[] valueCounter = new int[store.getNValues()];
		((KMerSortedArray<String>) store).visit(new KMerSortedArrayVisitor<String>() {
			@Override
			public void nextValue(KMerStore<String> trie, long kmer, short index, long i) {
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
}