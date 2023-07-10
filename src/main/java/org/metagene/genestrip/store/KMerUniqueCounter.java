package org.metagene.genestrip.store;

import it.unimi.dsi.fastutil.objects.Object2IntMap;

public interface KMerUniqueCounter {

	void clear();

	void put(long kmer, String taxid, long index);

	Object2IntMap<String> getUniqueKmerCounts();

}