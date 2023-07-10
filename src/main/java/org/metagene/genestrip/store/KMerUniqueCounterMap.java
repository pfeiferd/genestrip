package org.metagene.genestrip.store;

import it.unimi.dsi.fastutil.longs.Long2ObjectMap;
import it.unimi.dsi.fastutil.longs.Long2ObjectOpenHashMap;
import it.unimi.dsi.fastutil.objects.Object2IntLinkedOpenHashMap;
import it.unimi.dsi.fastutil.objects.Object2IntMap;

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
	public Object2IntMap<String> getUniqueKmerCounts() {
		Object2IntMap<String> res = new Object2IntLinkedOpenHashMap<String>();
		for (long kmer : map.keySet()) {
			String taxid = map.get(kmer);
			res.put(taxid, res.getInt(taxid) + 1);
		}

		return res;
	}
}