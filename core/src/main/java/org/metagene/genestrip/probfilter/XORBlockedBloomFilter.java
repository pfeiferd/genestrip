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
 * The {@code XOR} variant of {@link BlockedBloomFilter}: it replaces the MurmurHash3 finalizer of
 * {@link BlockedBloomFilter#hash(long)} by a bare exclusive or with the seed, and pays for that by
 * reducing with a modulo instead of the base class' multiply-shift.
 * <p>
 * The two go together. A multiply-shift is driven by the <em>high</em> bits of its input while a k-mer
 * carries its entropy in the low ones, so it is only sound after a hash that has carried that entropy
 * upwards. This hash does not, hence the modulo, which consumes every bit of the hash instead.
 * <p>
 * This filter suffers less from the mismatched pairing than {@link XORBloomFilter} would, because
 * it folds the hash by a 32 bit rotation for the bit positions anyway, but it does suffer: measured
 * over 1e6 k-mers at 10 bits per key, keeping the inherited multiply-shift costs 1.32% against this
 * variant's 1.30% at {@code k=31}, and 4.72% against 1.47% at {@code k=16}, where a key has only 32
 * significant bits. Random 64-bit keys hide the effect entirely, so measure on k-mer-like keys.
 * <p>
 * Whether this variant is also the faster one depends on the backing, because the two reduce in
 * different places. On the small {@code int}-indexed backing the mixing filter reduces by the cheap
 * 32-bit {@code reduceInt}, and the division of the modulo then outweighs what the simpler hash saves:
 * measured over 1e7 k-mers at 10 bits per key this variant looks up some 7 to 15% <em>slower</em> than
 * {@link BlockedBloomFilter}. On the bucketed backing the mixing filter reduces by
 * {@code Math.multiplyHigh} instead, which tips the balance: there this variant is about as fast on
 * hits and some 11 to 17% <em>faster</em> on misses. Unlike {@link XORBloomFilter}, which evaluates
 * its hash seven times per lookup, this filter hashes once, so the margins are small either way and
 * they were measured on aarch64 - re-measure before relying on them elsewhere.
 */
public class XORBlockedBloomFilter extends BlockedBloomFilter {
    private static final long serialVersionUID = 1L;

    /**
     * Creates a filter sized for {@code expectedInsertions} k-mers with
     * {@link #DEFAULT_BITS_PER_KEY} bits per key.
     *
     * @param expectedInsertions the expected number of k-mers to be inserted
     */
    public XORBlockedBloomFilter(long expectedInsertions) {
        super(expectedInsertions);
    }

    /**
     * Creates a filter sized for {@code expectedInsertions} k-mers with the given number of bits per
     * key and a fixed default seed.
     *
     * @param expectedInsertions the expected number of k-mers to be inserted
     * @param bitsPerKey the number of bits allocated per key
     */
    public XORBlockedBloomFilter(long expectedInsertions, int bitsPerKey) {
        super(expectedInsertions, bitsPerKey);
    }

    /**
     * Creates a filter sized for {@code expectedInsertions} k-mers with the given number of bits per
     * key and hash seed.
     *
     * @param expectedInsertions the expected number of k-mers to be inserted
     * @param bitsPerKey the number of bits allocated per key
     * @param seed       the hash seed used to derive bit positions
     */
    public XORBlockedBloomFilter(long expectedInsertions, int bitsPerKey, long seed) {
        super(expectedInsertions, bitsPerKey, seed);
    }

    /**
     * Creates a filter sized for {@code expectedInsertions} k-mers with the given number of bits per
     * key, hash seed and large-backing bucket width.
     *
     * @param expectedInsertions the expected number of k-mers to be inserted
     * @param bitsPerKey  the number of bits allocated per key
     * @param seed        the hash seed used to derive bit positions
     * @param bucketShift base-2 logarithm of the large-backing bucket width in words
     */
    public XORBlockedBloomFilter(long expectedInsertions, int bitsPerKey, long seed, int bucketShift) {
        super(expectedInsertions, bitsPerKey, seed, bucketShift);
    }

    /**
     * Creates a filter able to take the bucketed backing regardless of the filter's size, so that a test
     * can exercise that path at a size which fits in memory.
     *
     * @param expectedInsertions the expected number of k-mers to be inserted
     * @param bitsPerKey  the number of bits allocated per key
     * @param seed        the hash seed used to derive bit positions
     * @param bucketShift base-2 logarithm of the large-backing bucket width in words
     * @param forceLarge  whether to use the bucketed backing even when the small one would suffice
     */
    XORBlockedBloomFilter(long expectedInsertions, int bitsPerKey, long seed, int bucketShift,
                          boolean forceLarge) {
        super(expectedInsertions, bitsPerKey, seed, bucketShift, forceLarge);
    }

    /**
     * Creates a bucket-backed filter of this variant, for tests that need the bucketed backing at a size
     * that fits in memory.
     *
     * @param expectedInsertions the expected number of k-mers to be inserted
     * @param bitsPerKey the number of bits allocated per key
     * @return a filter of the given sizing that uses the bucketed backing
     */
    static XORBlockedBloomFilter newLargeBackedXOR(long expectedInsertions, int bitsPerKey) {
        return new XORBlockedBloomFilter(expectedInsertions, bitsPerKey, DEFAULT_SEED,
                minBucketShift(expectedInsertions, bitsPerKey), true);
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
     * {@link BlockedBloomFilter#reduce(long)}, because {@link #hash(long)} does not mix and a
     * multiply-shift would therefore build the word index from near-constant high bits.
     *
     * @param v the hash value to reduce
     * @return the start word index in {@code [0, buckets)} for the given hash value
     */
    @Override
    protected final long reduce(final long v) {
        return Math.abs(v % buckets);
    }

    /**
     * Reduces by a modulo, for the same reason {@link #reduce(long)} does. The cast is safe because the
     * small backing is only chosen while the word count fits an {@code int}.
     *
     * @param v the hash value to reduce
     * @return the start word index in {@code [0, buckets)} for the given hash value
     */
    @Override
    protected final int reduceInt(final long v) {
        return (int) Math.abs(v % buckets);
    }
}
