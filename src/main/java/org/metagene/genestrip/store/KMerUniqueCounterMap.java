package org.metagene.genestrip.store;

import it.unimi.dsi.fastutil.longs.Long2ObjectMap;
import it.unimi.dsi.fastutil.longs.Long2ObjectOpenHashMap;
import it.unimi.dsi.fastutil.objects.Object2LongLinkedOpenHashMap;
import it.unimi.dsi.fastutil.objects.Object2LongMap;

public class KMerUniqueCounterMap implements KMerUniqueCounter {
	private final Long2ObjectMap<String> map;

	public KMerUniqueCounterMap() {
		map = new Long2ObjectOpenHashMap<String>();
	}

	@Override	
	public void clear() {
		map.clear();
	}

	@Override	
	public synchronized void put(long kmer, String taxid, long index) {
		map.put(kmer, taxid);
	}

	@Override	
	public Object2LongMap<String> getUniqueKmerCounts() {
		Object2LongMap<String> res = new Object2LongLinkedOpenHashMap<String>();
		for (long kmer : map.keySet()) {
			String taxid = map.get(kmer);
			res.put(taxid, res.getLong(taxid) + 1);
		}

		return res;
	}
	
	@Override
	public int getUniqueKmerCount(String taxid) {
		int res = 0;
		for (long kmer : map.keySet()) {
			String taxid2 = map.get(kmer);
			if (taxid.equals(taxid2)) {
				res++;
			}
		}

		return res;
	}
}