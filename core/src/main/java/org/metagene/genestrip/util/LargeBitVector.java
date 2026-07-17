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

import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.Serializable;
import java.util.Arrays;

/**
 * A bit vector that transparently switches between a single {@code long[]} for small capacities and a
 * self-managed array of {@code long[]} buckets for capacities beyond the range of a normal array, and
 * never shrinks.
 * <p>
 * The large backing is a plain {@code long[][]} (an "array of arrays") laid out on a fixed grid: every
 * bucket holds exactly {@link #BUCKET_SIZE} words except the last, which holds the remainder. The grid
 * is self-managed (no external big-array library is involved); the bucket width matches the historic
 * segmentation, so the serialized form is unchanged. The buckets serve only to address storage beyond
 * a single {@code long[]}; they are made as few and large as the capacity requires and carry no lock.
 * <p>
 * Concurrent inserts go through {@link #setAndTestWasUnset(long)}, which is made thread-safe with a
 * fixed pool of {@link #SYNCS} stripe locks ({@link #syncs}) selected by the affected word index —
 * the same technique {@code AbstractKMerStore} uses for its entry updates. Decoupling the lock count
 * from the bucket count keeps concurrency high (up to {@link #SYNCS} writers in parallel) no matter
 * how few buckets the storage happens to use; two writers touching the same word always pick the same
 * stripe, so no update is lost.
 */
public class LargeBitVector implements Serializable {
    private static final long serialVersionUID = 1L;

    /** Maximum capacity in longs for which the small {@code long[]} backing is used. */
    public static long MAX_SMALL_CAPACITY = Integer.MAX_VALUE - 8;

    /** Base-2 logarithm of {@link #BUCKET_SIZE}; a word index is split into bucket and displacement here. */
    private static final int BUCKET_SHIFT = 27;
    /** Number of {@code long} words per bucket of the large backing (all but the last are full). */
    private static final int BUCKET_SIZE = 1 << BUCKET_SHIFT;
    /** Mask selecting the in-bucket displacement of a word index. */
    private static final int BUCKET_MASK = BUCKET_SIZE - 1;

    /** Number of stripe locks used to serialize concurrent word updates. Must be a power of two. */
    private static final int SYNCS = 512;
    /** Mask selecting a stripe lock from a word index. */
    private static final int SYNCS_MASK = SYNCS - 1;

    // All made public for inlining (optimization):
    /** The current capacity in 64-bit words. */
    public long size;
    /** The small backing array, or {@code null} when the large backing is used. */
    public long[] bits;
    /** The large (bucketed) backing, or {@code null} when the small backing is used. */
    public long[][] largeBits;
    /**
     * Stripe locks serializing concurrent {@link #setAndTestWasUnset(long)} calls, keyed by word
     * index. Transient (locks carry no state worth serializing) and rebuilt on construction and after
     * deserialization.
     */
    private transient Object[] syncs;

    /**
     * Creates a bit vector with at least the given initial capacity in bits.
     *
     * @param initialSize the initial capacity in bits
     */
    public LargeBitVector(long initialSize) {
        this(initialSize, false);
    }

    /**
     * Creates a bit vector with at least the given initial capacity in bits.
     *
     * @param initialSize  the initial capacity in bits
     * @param enforceLarge if true, forces the big-array backing even when a normal array would fit.
     */
    public LargeBitVector(long initialSize, boolean enforceLarge) {
        size = -1;
        initSyncs();
        ensureCapacity(initialSize, enforceLarge);
    }

    private void initSyncs() {
        syncs = new Object[SYNCS];
        for (int i = 0; i < syncs.length; i++) {
            syncs[i] = new Object();
        }
    }

    /**
     * Returns the stripe lock for the given word index, mapping it onto the fixed {@link #syncs} pool
     * with a cheap mask (the pool size is a power of two, so no modulo is needed). Writers to the same
     * word share a lock; writers to other words spread across the pool.
     *
     * @param wordIndex the affected word index
     * @return the stripe lock to synchronize on
     */
    private Object syncFor(long wordIndex) {
        return syncs[(int) (wordIndex & SYNCS_MASK)];
    }

    /**
     * Rebuilds the transient stripe locks after deserialization.
     *
     * @param ois the stream to read from
     * @throws IOException            if reading fails
     * @throws ClassNotFoundException if a serialized class cannot be resolved
     */
    private void readObject(ObjectInputStream ois) throws IOException, ClassNotFoundException {
        ois.defaultReadObject();
        initSyncs();
    }

    /**
     * Sets all bits to zero.
     */
    public void clear() {
        if (largeBits != null) {
            for (long[] bucket : largeBits) {
                Arrays.fill(bucket, 0L);
            }
        } else {
            Arrays.fill(bits, 0L);
        }
    }

    /**
     * Ensure that the bit vector has at least the desired size. The bit vector
     * might be enlarged but never gets smaller in size.
     *
     * @param newSize The desired size in bits.
     * @param enforceLarge if true, forces the big-array backing even when a normal array would fit.
     * @return Whether the vector got bigger or not.
     */
    public boolean ensureCapacity(long newSize, boolean enforceLarge) {
        newSize = (newSize + 63) / 64;
        if (newSize > size) {
            long oldWords = size < 0 ? 0 : size;
            size = newSize;
            if (size > MAX_SMALL_CAPACITY || enforceLarge || largeBits != null) {
                largeBits = growLarge(bits, largeBits, oldWords, size);
                bits = null;
            } else {
                long[] oldBits = bits;
                bits = new long[(int) size];
                if (oldBits != null) {
                    System.arraycopy(oldBits, 0, bits, 0, oldBits.length);
                }
            }
            return true;
        } else {
            return false;
        }
    }

    /**
     * Allocates the large (bucketed) backing for {@code newWords} words on the fixed {@link
     * #BUCKET_SIZE} grid and copies the first {@code oldWords} words of the previous backing into it.
     * Exactly one of {@code smallOld} / {@code largeOld} is non-null (the previous backing), or both
     * are null on the very first allocation.
     *
     * @param smallOld the previous small backing, or {@code null}
     * @param largeOld the previous large backing, or {@code null}
     * @param oldWords the number of words to preserve from the previous backing
     * @param newWords the desired capacity in words
     * @return the freshly allocated large backing holding the preserved words
     */
    private static long[][] growLarge(long[] smallOld, long[][] largeOld, long oldWords, long newWords) {
        int bucketCount = (int) ((newWords + BUCKET_SIZE - 1) >>> BUCKET_SHIFT);
        long[][] result = new long[bucketCount][];
        for (int b = 0; b < bucketCount; b++) {
            long start = (long) b << BUCKET_SHIFT;
            result[b] = new long[(int) Math.min(BUCKET_SIZE, newWords - start)];
        }
        // Copy the retained words bucket-aligned. Source and destination share the same grid, so a
        // block never straddles more than one source and one destination bucket.
        long copied = 0;
        while (copied < oldWords) {
            int dstBucket = (int) (copied >>> BUCKET_SHIFT);
            int dstOff = (int) (copied & BUCKET_MASK);
            long[] src;
            int srcOff;
            int srcRoom;
            if (largeOld != null) {
                src = largeOld[(int) (copied >>> BUCKET_SHIFT)];
                srcOff = (int) (copied & BUCKET_MASK);
                srcRoom = src.length - srcOff;
            } else {
                src = smallOld;
                srcOff = (int) copied;
                srcRoom = smallOld.length - srcOff;
            }
            int n = (int) Math.min(Math.min(result[dstBucket].length - dstOff, srcRoom), oldWords - copied);
            System.arraycopy(src, srcOff, result[dstBucket], dstOff, n);
            copied += n;
        }
        return result;
    }

    /**
     * Returns whether this vector currently uses the large (big-array) backing storage.
     *
     * @return {@code true} if the large backing is in use
     */
    public boolean isLarge() {
        return largeBits != null;
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
        if (largeBits != null) {
            final long arrayIndex = index >>> 6;
            largeBits[(int) (arrayIndex >>> BUCKET_SHIFT)][(int) (arrayIndex & BUCKET_MASK)] &= mask;
        } else {
            int arrayIndex = (int) (index >>> 6);
            bits[arrayIndex] = bits[arrayIndex] & mask;
        }
    }

    /**
     * Sets (to 1) the bit at the given index.
     *
     * @param index the bit index to set
     */
    public final void set(final long index) {
        if (largeBits != null) {
            final long arrayIndex = index >>> 6;
            largeBits[(int) (arrayIndex >>> BUCKET_SHIFT)][(int) (arrayIndex & BUCKET_MASK)] |= (1L << (index & 0b111111));
        } else {
            // Using optimization instead of '%', see:
            // http://lemire.me/blog/2016/06/27/a-fast-alternative-to-the-modulo-reduction/
            // and Line 34 in https://github.com/FastFilter/fastfilter_java/blob/master/fastfilter/src/main/java/org/fastfilter/utils/Hash.java
            // Not sure whether it would also work for long - probably not.
            //int arrayIndex = (int) ((((index >>> 6) & 0xffffffffL) * (size & 0xffffffffL)) >>> 32);
            // Original code:
            int arrayIndex = (int) (index >>> 6);
            bits[arrayIndex] = bits[arrayIndex] | (1L << (index & 0b111111));
        }
    }

    /**
     * Sets (to 1) the bit at the given index and reports whether this call was the one that set it,
     * i.e. whether the bit was previously 0. Unlike {@link #set(long)} this is safe for concurrent use
     * by multiple threads: the read-modify-write of the affected word runs under the stripe lock
     * selected by that word's index, so no concurrent bit set on the same word is ever lost while
     * writes to other words proceed in parallel (up to {@link #SYNCS} at a time).
     *
     * @param index the bit index to set
     * @return {@code true} if the bit transitioned from 0 to 1, {@code false} if it was already set
     */
    public final boolean setAndTestWasUnset(final long index) {
        final long bit = 1L << (index & 0b111111);
        final long arrayIndex = index >>> 6;
        final Object lock = syncFor(arrayIndex);
        if (largeBits != null) {
            final long[] bucket = largeBits[(int) (arrayIndex >>> BUCKET_SHIFT)];
            final int displacement = (int) (arrayIndex & BUCKET_MASK);
            synchronized (lock) {
                final long old = bucket[displacement];
                bucket[displacement] = old | bit;
                return (old & bit) == 0;
            }
        } else {
            final int i = (int) arrayIndex;
            synchronized (lock) {
                final long old = bits[i];
                bits[i] = old | bit;
                return (old & bit) == 0;
            }
        }
    }

    /**
     * Returns whether the bit at the given index is set.
     *
     * @param index the bit index to test
     * @return {@code true} if the bit is set
     */
    public final boolean get(final long index) {
        if (largeBits != null) {
            final long arrayIndex = index >>> 6;
            final long word = largeBits[(int) (arrayIndex >>> BUCKET_SHIFT)][(int) (arrayIndex & BUCKET_MASK)];
            return ((word >> (index & 0b111111)) & 1L) == 1;
        } else {
            // Using optimization instead of '%', see:
            // http://lemire.me/blog/2016/06/27/a-fast-alternative-to-the-modulo-reduction/
            // and Line 34 in https://github.com/FastFilter/fastfilter_java/blob/master/fastfilter/src/main/java/org/fastfilter/utils/Hash.java
            // Not sure whether it would also work for long - probably not.
            //int arrayIndex = (int) ((((index >>> 6) & 0xffffffffL) * (size & 0xffffffffL)) >>> 32);
            // Original code:
            int arrayIndex = (int) (index >>> 6);
            return (((bits[arrayIndex] >> (index & 0b111111)) & 1L) == 1);
        }
    }
}
