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