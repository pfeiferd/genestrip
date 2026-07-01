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

import it.unimi.dsi.fastutil.objects.Object2LongMap;

/**
 * Tracks distinct k-mers matched (looked up) per taxid, so the number of unique matched k-mers can be
 * reported per taxon.
 */
public interface KMerUniqueCounter {
	/**
	 * Clears all recorded unique-k-mer counts.
	 */
	public void clear();

	/**
	 * Records that the given k-mer, stored at position {@code index} in the backing store, was matched
	 * for the given taxid.
	 *
	 * @param kmer the matched k-mer
	 * @param taxid the taxid the k-mer was matched for
	 * @param index the position of the k-mer in the backing store
	 */
	public void put(long kmer, String taxid, long index);

	/**
	 * Returns the number of distinct matched k-mers per taxid.
	 *
	 * @return the number of distinct matched k-mers per taxid.
	 */
	public Object2LongMap<String> getUniqueKmerCounts();

	/**
	 * Returns the number of distinct matched k-mers for the given taxid.
	 *
	 * @param taxid the taxid to query
	 * @return the number of distinct matched k-mers for the given taxid.
	 */
	public int getUniqueKmerCount(String taxid);
}