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
package org.metagene.genestrip.refseq;

import org.metagene.genestrip.tax.TaxTree.TaxIdNode;

/**
 * Maps genomic/RNA accession numbers to the tax id nodes they belong to.
 */
public interface AccessionMap {
	/**
	 * Adds a mapping from the accession key {@code array[start, end)} to the given tax id node.
	 *
	 * @param array the byte array containing the accession key.
	 * @param start the start index of the key (inclusive).
	 * @param end the end index of the key (exclusive).
	 * @param node the tax id node to map the accession key to.
	 */
	public void put(byte[] array, int start, int end, TaxIdNode node);
	/**
	 * Returns the tax id node for the accession key {@code array[start, end)}, or null if none.
	 *
	 * @param array the byte array containing the accession key.
	 * @param start the start index of the key (inclusive).
	 * @param end the end index of the key (exclusive).
	 * @param completeGenomesOnly if true, only complete-genome, RNA and mRNA accessions are resolved.
	 * @return the tax id node for the accession key, or null if none.
	 */
	public TaxIdNode get(byte[] array, int start, int end, boolean completeGenomesOnly);
	/**
	 * Prepares the map for lookups; must be called after all {@link #put} calls and before any
	 * {@link #get}.
	 */
	public void optimize();
	/**
	 * Returns the number of accession entries mapped to the given node.
	 *
	 * @param node the tax id node to count entries for.
	 * @return the number of accession entries mapped to the node.
	 */
	public int getEntriesForNode(TaxIdNode node);
}
