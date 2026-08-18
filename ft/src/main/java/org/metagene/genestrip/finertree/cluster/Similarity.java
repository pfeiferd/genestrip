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
package org.metagene.genestrip.finertree.cluster;

/**
 * Provides pairwise similarity values over a fixed set of items, indexed from {@code 0} to
 * {@link #values()}{@code  - 1}. Used as the input to hierarchical agglomerative clustering in
 * {@link SimpleAggloClustering}.
 */
public interface Similarity {
    /**
     * Returns the number of items over which similarities are defined.
     *
     * @return the item count
     */
    public int values();

    /**
     * Returns the similarity between the items at indices {@code i} and {@code j}. Implementations
     * are expected to be symmetric, i.e. {@code getSimilarity(i, j) == getSimilarity(j, i)}.
     *
     * @param i the index of the first item
     * @param j the index of the second item
     * @return the similarity between the two items
     */
    public double getSimilarity(int i, int j);
}
