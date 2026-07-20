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
 * A bit vector backed by a self-managed array of {@code long[]} buckets that never shrinks.
 * <p>
 * The backing is a plain {@code long[][]} (an "array of arrays") laid out on a fixed, power-of-two
 * grid: every bucket holds exactly {@code 1 << bucketShift} words except the last, which holds the
 * remainder. The grid is self-managed (no external big-array library is involved). The bucket width is
 * chosen once, at the first allocation, so that the vector holds at least {@code 2^minBucketBits}
 * buckets (capped at the vector's word count) while keeping the buckets as few and large as that lower
 * bound allows, up to {@link #MAX_BUCKET_SHIFT} words each.
 * <p>
 * Concurrent inserts go through {@link #setAndTestWasUnset(long)}, which is made thread-safe by
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
     * Default {@link #minBucketBits}: {@code 9}, i.e. a floor of {@code 2^9 == 512} buckets — enough
     * independent locks for efficient multi-threaded inserts.
     */
    public static final int DEFAULT_MIN_BUCKET_BITS = 9;

    /** Upper bound on {@link #minBucketBits}, keeping {@code 1 << minBucketBits} a sane bucket floor. */
    private static final int MAX_MIN_BUCKET_BITS = 30;

    /**
     * Maximum base-2 logarithm of the bucket width, i.e. the largest bucket holds {@code 1 << 27}
     * words. Bounds each bucket's length to a comfortably int-addressable size and caps the outer
     * array's growth; it also matches the historic bucket width.
     */
    private static final int MAX_BUCKET_SHIFT = 27;

    /** The current capacity in 64-bit words. */
    private long size;
    /**
     * Base-2 logarithm of the lower bound on the number of buckets: the vector holds at least
     * {@code 1 << minBucketBits} buckets (hence independent locks), capped at its word count. Expressing
     * the floor as a bit count keeps {@link #computeBucketShift} a plain subtraction of bit widths.
     */
    private final int minBucketBits;
    /**
     * Base-2 logarithm of the bucket width; a word index is split into bucket ({@code >>> bucketShift})
     * and in-bucket displacement ({@code & bucketMask}). {@code -1} until the first non-empty
     * allocation fixes it; unchanged afterwards so the grid stays stable as the vector grows.
     */
    private int bucketShift;
    /** Mask selecting the in-bucket displacement of a word index ({@code (1 << bucketShift) - 1}). */
    private int bucketMask;
    /** The bucketed backing. */
    private long[][] largeBits;
    /**
     * Per-bucket count of bits set through {@link #setAndTestWasUnset(long)}, one entry per bucket (so
     * {@code setBitsPerBucket.length == largeBits.length}). Each counter is only ever touched under its
     * own bucket's monitor, so the concurrent set path updates it without contending across buckets;
     * {@link #getBitsEverSet()} sums the array. See that method for what is and isn't counted.
     */
    private long[] setBitsPerBucket;

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
        size = -1;
        this.minBucketBits = minBucketBits;
        bucketShift = -1;
        ensureCapacity(initialSize);
    }

    /**
     * Sets all bits to zero.
     */
    public void clear() {
        for (long[] bucket : largeBits) {
            Arrays.fill(bucket, 0L);
        }
        Arrays.fill(setBitsPerBucket, 0L);
    }

    /**
     * Ensure that the bit vector has at least the desired size. The bit vector
     * might be enlarged but never gets smaller in size.
     *
     * @param newSize The desired size in bits.
     * @return Whether the vector got bigger or not.
     */
    public boolean ensureCapacity(long newSize) {
        newSize = (newSize + 63) / 64;
        if (newSize > size) {
            long oldWords = size < 0 ? 0 : size;
            size = newSize;
            // Fix the bucket grid on the first non-empty allocation and keep it, so growth only adds
            // buckets on the same grid (making the retained words a straight bucket-aligned copy).
            if (bucketShift < 0 && newSize > 0) {
                bucketShift = computeBucketShift(newSize, minBucketBits);
                bucketMask = (1 << bucketShift) - 1;
            }
            largeBits = growLarge(largeBits, oldWords, size);
            // Grow the per-bucket counters to match; growth keeps the grid, so bucket b still holds the
            // same preserved words - its count carries over, and the freshly added buckets start at 0.
            long[] newCounters = new long[largeBits.length];
            if (setBitsPerBucket != null) {
                System.arraycopy(setBitsPerBucket, 0, newCounters, 0, setBitsPerBucket.length);
            }
            setBitsPerBucket = newCounters;
            return true;
        } else {
            return false;
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
     * Allocates the bucketed backing for {@code newWords} words on the fixed grid ({@link #bucketShift})
     * and copies the first {@code oldWords} words of the previous backing into it. {@code largeOld} is
     * the previous backing, or {@code null} on the very first allocation.
     *
     * @param largeOld the previous backing, or {@code null}
     * @param oldWords the number of words to preserve from the previous backing
     * @param newWords the desired capacity in words
     * @return the freshly allocated backing holding the preserved words
     */
    private long[][] growLarge(long[][] largeOld, long oldWords, long newWords) {
        if (newWords == 0) {
            return new long[0][];
        }
        final int bucketSize = 1 << bucketShift;
        int bucketCount = (int) bucketCount(newWords, bucketShift);
        long[][] result = new long[bucketCount][];
        for (int b = 0; b < bucketCount; b++) {
            long start = (long) b << bucketShift;
            result[b] = new long[(int) Math.min(bucketSize, newWords - start)];
        }
        // Copy the retained words bucket-aligned. The old and new backings share the same grid (the
        // shift is fixed after the first allocation), so a block never straddles more than one source
        // and one destination bucket.
        long copied = 0;
        while (copied < oldWords) {
            int bucket = (int) (copied >>> bucketShift);
            int off = (int) (copied & bucketMask);
            long[] src = largeOld[bucket];
            int srcRoom = src.length - off;
            int n = (int) Math.min(Math.min(result[bucket].length - off, srcRoom), oldWords - copied);
            System.arraycopy(src, off, result[bucket], off, n);
            copied += n;
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
     * Clears (sets to 0) the bit at the given index.
     *
     * @param index the bit index to clear
     */
    public void clear(long index) {
        final long mask = ~(1L << (index & 0b111111));
        final long arrayIndex = index >>> 6;
        largeBits[(int) (arrayIndex >>> bucketShift)][(int) (arrayIndex & bucketMask)] &= mask;
    }

    /**
     * Sets (to 1) the bit at the given index.
     *
     * @param index the bit index to set
     */
    public final void set(final long index) {
        final long arrayIndex = index >>> 6;
        final int bucketIndex = (int) (arrayIndex >>> bucketShift);
        largeBits[bucketIndex][(int) (arrayIndex & bucketMask)] |= (1L << (index & 0b111111));
        setBitsPerBucket[bucketIndex]++;
    }

    /**
     * Sets (to 1) the bit at the given index and reports whether this call was the one that set it,
     * i.e. whether the bit was previously 0. Unlike {@link #set(long)} this is safe for concurrent use
     * by multiple threads: the read-modify-write of the affected word runs under the lock of the bucket
     * that owns it, so no concurrent bit set on the same word is ever lost while writes to other
     * buckets proceed in parallel.
     *
     * @param index the bit index to set
     * @return {@code true} if the bit transitioned from 0 to 1, {@code false} if it was already set
     */
    public final boolean setAndTestWasUnset(final long index) {
        final long bit = 1L << (index & 0b111111);
        final long arrayIndex = index >>> 6;
        final int bucketIndex = (int) (arrayIndex >>> bucketShift);
        final long[] bucket = largeBits[bucketIndex];
        final int displacement = (int) (arrayIndex & bucketMask);
        synchronized (bucket) {
            // Count every set (like set() does), still under the bucket's monitor so the shared
            // counter never races; this makes getBitsEverSet() total the number of set operations.
            setBitsPerBucket[bucketIndex]++;
            final long old = bucket[displacement];
            bucket[displacement] = old | bit;
            return (old & bit) == 0;
        }
    }

    /**
     * Returns the number of set operations performed, i.e. every {@link #set(long)} and
     * {@link #setAndTestWasUnset(long)} call is counted (even one that re-sets an already-set bit);
     * cleared bits are not subtracted. So this is the total count of bit-set operations, not the number
     * of distinct bits currently set.
     *
     * @return the number of set operations performed
     */
    public long getBitsEverSet() {
        long total = 0;
        for (long count : setBitsPerBucket) {
            total += count;
        }
        return total;
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
