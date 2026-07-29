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
     * Multiplier spreading the index over the whole key. Any odd constant does, as multiplying by one
     * is a bijection modulo 2^64 and hence maps distinct indexes to distinct masks. It deliberately
     * differs from the constant {@link org.metagene.genestrip.bloom.BlockedKMerBloomFilter} multiplies
     * by when deriving a word position, so that the two steps cannot interact.
     */
    private static final long INDEX_MULTIPLIER = 0xD6E8FEB86659FD93L;

    /**
     * Folds a child index into a k-mer so that the resulting pair can be stored in and queried from
     * a filter over {@code long} keys. The index is spread over all 64 bits by
     * {@link #INDEX_MULTIPLIER} before being xored in, so that a k-mer combined with different
     * indexes yields keys that differ in every part of the word.
     * <p>
     * That the whole word is affected is what the filters this feeds depend on: they hash a key
     * trivially as {@code seed ^ key} and therefore rely on the key itself to carry the entropy.
     * Xoring the index in unmixed - once into each half of the key - yields a mask that is invariant
     * under a 32 bit rotation. Such a mask cancels out entirely where the filter folds the key by that
     * rotation to derive the bit positions, and its low bits are dropped again where the filter derives
     * the word position from the upper half of a product. The index would then hardly reach the filter
     * at all, and k-mers stored under one index would be reported as present under any other.
     *
     * @param data  the k-mer, encoded as a {@code long}
     * @param index the index of the child subtree the k-mer is attributed to
     * @return the combined key identifying the (k-mer, child index) pair
     */
    public static long combine(final long data, final int index) {
        return data ^ (index * INDEX_MULTIPLIER);
    }
}
