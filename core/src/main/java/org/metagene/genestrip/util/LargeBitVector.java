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
package org.metagene.genestrip.util;

import java.io.Serializable;
import java.util.Arrays;

/**
 * A bit vector backed by a self-managed array of {@code long[]} buckets whose size is fixed at
 * construction.
 * <p>
 * The backing is a plain {@code long[][]} (an "array of arrays") laid out on a fixed, power-of-two
 * grid: every bucket holds exactly {@code 1 << bucketShift} words except the last, which holds the
 * remainder. The grid is self-managed (no external big-array library is involved). The bucket width is
 * chosen once, at construction, so that the vector holds at least {@code 2^minBucketBits}
 * buckets (capped at the vector's word count) while keeping the buckets as few and large as that lower
 * bound allows, up to {@link #MAX_BUCKET_SHIFT} words each.
 * <p>
 * Concurrent inserts go through {@link #set(long)}, which is made thread-safe by
 * synchronizing on the bucket that owns the affected word — the same technique {@code RadixKMerStore}
 * uses for its entry updates. The bucket is both the storage and the lock, so no separate stripe key
 * has to be computed: two writers touching the same word share the same bucket (hence the same lock),
 * while writers to different buckets proceed in parallel. The {@code 2^minBucketBits} floor therefore
 * doubles as the lower bound on the number of independent locks, keeping concurrency high even for
 * modestly sized vectors; the default floor is {@link #DEFAULT_MIN_BUCKET_BITS} bits.
 */
public class LargeBitVector implements Serializable {
    private static final long serialVersionUID = 2L;

    /**
     * Default buckets: {@code 9}, i.e. a floor of {@code 2^9 == 512} buckets — enough
     * independent locks for efficient multi-threaded inserts.
     */
    public static final int DEFAULT_MIN_BUCKET_BITS = 9;

    private static final int MAX_MIN_BUCKET_BITS = 30;

    /**
     * Maximum base-2 logarithm of the bucket width, i.e. the largest bucket holds {@code 1 << 27}
     * words. Bounds each bucket's length to a comfortably int-addressable size and caps the outer
     * array's growth; it also matches the historic bucket width.
     */
    private static final int MAX_BUCKET_SHIFT = 27;

    /** The capacity in 64-bit words, fixed at construction. */
    private final long size;

    /**
     * Base-2 logarithm of the bucket width; a word index is split into bucket ({@code >>> bucketShift})
     * and in-bucket displacement ({@code & bucketMask}). {@code -1} for an empty vector; otherwise chosen
     * at construction and fixed for the vector's lifetime.
     */
    private final int bucketShift;
    /** Mask selecting the in-bucket displacement of a word index ({@code (1 << bucketShift) - 1}). */
    private final int bucketMask;
    /** The bucketed backing, allocated once at construction. */
    private final long[][] largeBits;

    /**
     * Creates a bit vector with at least the given initial capacity in bits and the default bucket
     * floor ({@link #DEFAULT_MIN_BUCKET_BITS}).
     *
     * @param initialSize the initial capacity in bits
     */
    public LargeBitVector(long initialSize) {
        this(initialSize, DEFAULT_MIN_BUCKET_BITS);
    }

    /**
     * Creates a bit vector with at least the given initial capacity in bits and at least
     * {@code 1 << minBucketBits} buckets (once non-empty), so that concurrent inserts have at least that
     * many independent locks to spread across.
     *
     * @param initialSize   the initial capacity in bits
     * @param minBucketBits base-2 logarithm of the minimum bucket (and lock) count; in
     *                      {@code [1, 30]}
     */
    public LargeBitVector(long initialSize, int minBucketBits) {
        if (minBucketBits < 1 || minBucketBits > MAX_MIN_BUCKET_BITS) {
            throw new IllegalArgumentException(
                    "minBucketBits must be in [1, " + MAX_MIN_BUCKET_BITS + "], got " + minBucketBits);
        }
        size = (initialSize + 63) / 64;
        // Fix the bucket grid from the (non-empty) capacity; an empty vector needs no grid.
        if (size > 0) {
            bucketShift = computeBucketShift(size, minBucketBits);
            bucketMask = (1 << bucketShift) - 1;
        } else {
            bucketShift = -1;
            bucketMask = 0;
        }
        largeBits = allocateBuckets(size);
    }

    /**
     * Sets all bits to zero.
     */
    public void clear() {
        for (long[] bucket : largeBits) {
            Arrays.fill(bucket, 0L);
        }
    }

    /**
     * Chooses the bucket-width exponent for a vector of {@code words} words: the largest exponent (up
     * to {@link #MAX_BUCKET_SHIFT}) that still yields at least {@code 2^minBucketBits} buckets, so the
     * buckets are as few and large as the lock floor permits. Because both the floor and the bucket
     * width are powers of two, this is just {@code floor(log2(words)) - minBucketBits}, clamped: a
     * vector of {@code < 2^minBucketBits} words gets one word per bucket (the most locks it can support).
     *
     * @param words         the capacity in words (must be positive)
     * @param minBucketBits base-2 logarithm of the desired minimum bucket count
     * @return the bucket-width exponent to use
     */
    private static int computeBucketShift(long words, int minBucketBits) {
        int wordsLog2 = 63 - Long.numberOfLeadingZeros(words); // floor(log2(words)), words > 0
        int shift = wordsLog2 - minBucketBits;
        if (shift < 0) {
            return 0;
        }
        return Math.min(shift, MAX_BUCKET_SHIFT);
    }

    /** Number of buckets of width {@code 1 << shift} needed to hold {@code words} words. */
    private static long bucketCount(long words, int shift) {
        return (words + (1L << shift) - 1) >>> shift;
    }

    /**
     * Allocates the zeroed bucketed backing for {@code words} words on the fixed grid
     * ({@link #bucketShift}).
     *
     * @param words the desired capacity in words
     * @return the freshly allocated backing
     */
    private long[][] allocateBuckets(long words) {
        if (words == 0) {
            return new long[0][];
        }
        final int bucketSize = 1 << bucketShift;
        int bucketCount = (int) bucketCount(words, bucketShift);
        long[][] result = new long[bucketCount][];
        for (int b = 0; b < bucketCount; b++) {
            long start = (long) b << bucketShift;
            result[b] = new long[(int) Math.min(bucketSize, words - start)];
        }
        return result;
    }

    /**
     * Returns the current capacity in bits.
     *
     * @return the current capacity in bits
     */
    public long getBitSize() {
        return size * 64;
    }

    /**
     * Clears (sets to 0) the bit at the given index and reports whether this call was the one that
     * cleared it, i.e. whether the bit was previously 1. Safe for concurrent use by multiple threads:
     * the read-modify-write of the affected word runs under the lock of the bucket that owns it, so no
     * concurrent update to the same word is ever lost while writes to other buckets proceed in parallel.
     *
     * @param index the bit index to clear
     * @return {@code true} if the bit transitioned from 1 to 0, {@code false} if it was already clear
     */
    public final boolean clear(final long index) {
        final long bit = 1L << (index & 0b111111);
        final long arrayIndex = index >>> 6;
        final int bucketIndex = (int) (arrayIndex >>> bucketShift);
        final long[] bucket = largeBits[bucketIndex];
        final int displacement = (int) (arrayIndex & bucketMask);
        synchronized (bucket) {
            final long old = bucket[displacement];
            bucket[displacement] = old & ~bit;
            return (old & bit) != 0;
        }
    }

    /**
     * Sets (to 1) the bit at the given index and reports whether this call was the one that set it,
     * i.e. whether the bit was previously 0. Safe for concurrent use by multiple threads: the
     * read-modify-write of the affected word runs under the lock of the bucket that owns it, so no
     * concurrent bit set on the same word is ever lost while writes to other buckets proceed in parallel.
     *
     * @param index the bit index to set
     * @return {@code true} if the bit transitioned from 0 to 1, {@code false} if it was already set
     */
    public final boolean set(final long index) {
        final long bit = 1L << (index & 0b111111);
        final long arrayIndex = index >>> 6;
        final int bucketIndex = (int) (arrayIndex >>> bucketShift);
        final long[] bucket = largeBits[bucketIndex];
        final int displacement = (int) (arrayIndex & bucketMask);
        synchronized (bucket) {
            final long old = bucket[displacement];
            bucket[displacement] = old | bit;
            return (old & bit) == 0;
        }
    }

    /**
     * Returns whether the bit at the given index is set.
     *
     * @param index the bit index to test
     * @return {@code true} if the bit is set
     */
    public final boolean get(final long index) {
        final long arrayIndex = index >>> 6;
        final long word = largeBits[(int) (arrayIndex >>> bucketShift)][(int) (arrayIndex & bucketMask)];
        return ((word >> (index & 0b111111)) & 1L) == 1;
    }
}
