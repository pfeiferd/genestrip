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

import org.metagene.genestrip.probfilter.BlockedBloomFilter;

/**
 * Static helper for the k-mer index Bloom filter of the finer-tree (FT) extension, which records
 * (k-mer, child index) pairs rather than plain k-mers.
 */
public class KMerIndexFilterHelper {
    // Static utility class - not meant to be instantiated.
    private KMerIndexFilterHelper() {
    }

    /**
     * Multiplier spreading the index over the whole key. Any odd constant does, as multiplication by an
     * odd number is a bijection modulo 2^64 and hence maps distinct indexes to distinct masks. It is
     * unrelated to what {@link BlockedBloomFilter} multiplies by when deriving a word position - that is
     * the filter's own word count, not a constant - so the two steps cannot interact.
     */
    private static final long INDEX_MULTIPLIER = 0xD6E8FEB86659FD93L;

    /**
     * Folds a child index into a k-mer so that the resulting pair can be stored in and queried from
     * a filter over {@code long} keys. The index is spread over all 64 bits by
     * {@link #INDEX_MULTIPLIER} before being xored in, so that a k-mer combined with different
     * indexes yields keys that differ in every part of the word.
     * <p>
     * No filter of the {@code probfilter} package depends on that today: the mixing ones run a
     * MurmurHash3 finalizer over the key, and the {@code XOR} ones, which do hash it trivially as
     * {@code seed ^ key}, reduce by a modulo that consumes every bit of the result. Either way the index
     * reaches the filter however it is folded in, and xoring it in unmixed measures the same
     * false-positive rate.
     * <p>
     * The multiplication is kept because it costs one instruction and removes the assumption. A filter
     * pairing a non-mixing hash with a multiply-shift reduction - the pairing {@code HashReducePairingTest}
     * forbids, but which a new implementation could reintroduce - takes a key's bit positions from the key
     * folded by a 32 bit rotation and its word position from the upper half of a product. An index xored
     * in once into each half of the key is invariant under that rotation and cancels in the fold, and one
     * left in the low bits is dropped by the product; k-mers stored under one index then get reported as
     * present under every other, at 86% cross index false positives against 1.3% with the multiplication.
     * <p>
     * A plain {@code data ^ index} fails for a second and simpler reason: an index bounded by
     * {@link org.metagene.genestrip.store.KMerSortedArray#MAX_VALUES} occupies the low 16 bits only and
     * never reaches the upper half of the key at all. Placing it there by a shift instead would collide
     * as soon as a store offers more values than a shift of that width can hold.
     *
     * @param data  the k-mer, encoded as a {@code long}
     * @param index the index of the child subtree the k-mer is attributed to
     * @return the combined key identifying the (k-mer, child index) pair
     */
    public static long combine(final long data, final int index) {
        return data ^ (index * INDEX_MULTIPLIER);
    }
}
