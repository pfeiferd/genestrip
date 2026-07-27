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
package org.metagene.genestrip.finertree.bloom;

/**
 * Static helper for the k-mer index Bloom filter of the finer-tree (FT) extension, which records
 * (k-mer, child index) pairs rather than plain k-mers.
 */
public class KMerIndexFilterHelper {
    // Static utility class - not meant to be instantiated.
    private KMerIndexFilterHelper() {
    }

    /**
     * Folds a child index into a k-mer so that the resulting pair can be stored in and queried from
     * a filter over {@code long} keys. The index is mixed into both halves of the key, so that a
     * k-mer combined with different indexes yields different keys across the full 64 bit range. The
     * mapping is its own inverse with respect to the index, i.e. combining twice with the same index
     * restores the original k-mer.
     *
     * @param data  the k-mer, encoded as a {@code long}
     * @param index the index of the child subtree the k-mer is attributed to
     * @return the combined key identifying the (k-mer, child index) pair
     */
    public static long combine(final long data, final int index) {
        return data ^ ((long) index) ^ (((long) index) << 32);
    }
}
