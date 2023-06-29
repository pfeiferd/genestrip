package org.metagene.genestrip.trie;

import java.io.Serializable;

public interface KMerStoreFactory {
	public <V extends Serializable> KMerStore<V> createKMerStore(Class<V> clazz, Object... params);
}
