package org.metagene.genestrip.store;

import it.unimi.dsi.fastutil.objects.Object2LongMap;

public interface KMerUniqueCounter {
	public void clear();

	public void put(long kmer, String taxid, long index);

	public Object2LongMap<String> getUniqueKmerCounts();
	
	public int getUniqueKmerCount(String taxid);
}