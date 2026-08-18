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
package org.metagene.genestrip.probfilter;

/**
 * The {@code XOR} variant of {@link SingleWordBloomFilter}: it replaces the MurmurHash3 finalizer
 * of {@link SingleWordBloomFilter#hash(long)} by a bare exclusive or with the seed, and pays for
 * that by reducing with a modulo instead of the base class' multiply-shift.
 * <p>
 * The two go together. A multiply-shift is driven by the <em>high</em> bits of its input while a k-mer
 * carries its entropy in the low ones, so it is only sound after a hash that has carried that entropy
 * upwards. This hash does not, hence the modulo, which consumes every bit of the hash instead.
 * <p>
 * As in {@link XORBlockedBloomFilter} the mismatched pairing hurts less here than in
 * {@link XORBloomFilter}, because the bit positions are folded out of the hash anyway, but it does
 * hurt: measured over 1e6 k-mers at 10 bits per key, keeping the inherited multiply-shift costs 1.79%
 * against this variant's 1.80% at {@code k=31}, and 5.33% against 1.92% at {@code k=16}, where a key
 * has only 32 significant bits. Random 64-bit keys hide the effect entirely, so measure on k-mer-like
 * keys.
 * <p>
 * As in {@link XORBlockedBloomFilter} the speed of this variant depends on the backing, since that
 * decides which reduction the mixing filter it is compared against uses. Measured over 1e7 k-mers at 10
 * bits per key, this variant looks up some 2 to 9% <em>slower</em> on the small {@code int}-indexed
 * backing, whose {@code reduceInt} is cheap, and some 11 to 21% <em>faster</em> on the bucketed one,
 * where the mixing filter pays a {@code Math.multiplyHigh}. The margins were measured on aarch64 -
 * re-measure before relying on them elsewhere.
 */
public class XORSingleWordBloomFilter extends SingleWordBloomFilter {
    private static final long serialVersionUID = 1L;

    /**
     * Creates a filter sized for {@code expectedInsertions} k-mers with
     * {@link #DEFAULT_BITS_PER_KEY} bits per key.
     *
     * @param expectedInsertions the expected number of k-mers to be inserted
     */
    public XORSingleWordBloomFilter(long expectedInsertions) {
        super(expectedInsertions);
    }

    /**
     * Creates a filter sized for {@code expectedInsertions} k-mers with the given number of bits per
     * key, setting {@link #optimalHashBits(int)} bits per key.
     *
     * @param expectedInsertions the expected number of k-mers to be inserted
     * @param bitsPerKey         the number of bits allocated per key
     */
    public XORSingleWordBloomFilter(long expectedInsertions, int bitsPerKey) {
        super(expectedInsertions, bitsPerKey);
    }

    /**
     * Creates a filter sized for {@code expectedInsertions} k-mers that sets the given number of bits
     * per key, with a fixed default seed.
     *
     * @param expectedInsertions the expected number of k-mers to be inserted
     * @param bitsPerKey         the number of bits allocated per key
     * @param hashBits           the number of bits a key sets within its word
     */
    public XORSingleWordBloomFilter(long expectedInsertions, int bitsPerKey, int hashBits) {
        super(expectedInsertions, bitsPerKey, hashBits);
    }

    /**
     * Creates a filter with the given sizing, seed and large-backing bucket width, able to take the
     * bucketed backing regardless of the filter's size.
     *
     * @param expectedInsertions the expected number of k-mers to be inserted
     * @param bitsPerKey         the number of bits allocated per key
     * @param hashBits           the number of bits a key sets within its word
     * @param seed               the hash seed
     * @param bucketShift        base-2 logarithm of the large-backing bucket width in words
     * @param forceLarge         whether to use the bucketed backing even when the small one would suffice
     */
    XORSingleWordBloomFilter(long expectedInsertions, int bitsPerKey, int hashBits, long seed, int bucketShift,
                             boolean forceLarge) {
        super(expectedInsertions, bitsPerKey, hashBits, seed, bucketShift, forceLarge);
    }

    /**
     * Creates a bucket-backed filter of this variant, for tests that need the bucketed backing at a size
     * that fits in memory.
     *
     * @param expectedInsertions the expected number of k-mers to be inserted
     * @param bitsPerKey the number of bits allocated per key
     * @return a filter of the given sizing that uses the bucketed backing
     */
    static XORSingleWordBloomFilter newLargeBackedXOR(long expectedInsertions, int bitsPerKey) {
        return new XORSingleWordBloomFilter(expectedInsertions, bitsPerKey, optimalHashBits(bitsPerKey),
                DEFAULT_SEED, minBucketShift(expectedInsertions, bitsPerKey), true);
    }

    /**
     * XORs the key with the seed - no mixing whatsoever, which is what makes it fast.
     *
     * @param x the key to hash
     * @return the hash of the given key
     */
    @Override
    protected final long hash(long x) {
        return seed ^ x;
    }

    /**
     * Reduces by a modulo rather than by the multiply-shift of
     * {@link SingleWordBloomFilter#reduce(long)}, because {@link #hash(long)} does not mix and a
     * multiply-shift would therefore build the word index from near-constant high bits.
     *
     * @param v the hash value to reduce
     * @return the word index in {@code [0, words)} for the given hash value
     */
    @Override
    protected final long reduce(final long v) {
        return Math.abs(v % words);
    }

    /**
     * Reduces by a modulo, for the same reason {@link #reduce(long)} does. The cast is safe because the
     * small backing is only chosen while the word count fits an {@code int}.
     *
     * @param v the hash value to reduce
     * @return the word index in {@code [0, words)} for the given hash value
     */
    @Override
    protected final int reduceInt(final long v) {
        return (int) Math.abs(v % words);
    }
}
