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
package org.metagene.genestrip.bloom;

/**
 * A {@link KMerProbFilter} that additionally tracks the number of k-mers added to it. Not every
 * filter counts its entries (e.g. {@link BlockedKMerBloomFilter} does not), so this capability is
 * split out from the base {@link KMerProbFilter} interface.
 */
public interface CountingKMerProbFilter extends KMerProbFilter {
    /**
     * Returns the number of k-mers that have been added to the filter.
     * <p>
     * For the concurrent implementations, two threads inserting the same absent k-mer may both
     * observe it as new, so under concurrent use this count may marginally over-count; membership
     * answers are never affected.
     *
     * @return the number of inserted entries
     */
    public long getEntries();
}
