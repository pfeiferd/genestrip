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
import java.lang.invoke.MethodHandles;
import java.lang.invoke.VarHandle;
import java.util.Arrays;

import it.unimi.dsi.fastutil.BigArrays;
import it.unimi.dsi.fastutil.longs.LongBigArrays;

/**
 * A bit vector that transparently switches between a single {@code long[]} for small capacities and
 * fastutil big arrays for capacities beyond the range of a normal array, and never shrinks.
 */
public class LargeBitVector implements Serializable {
    private static final long serialVersionUID = 1L;

    /** Maximum capacity in longs for which the small {@code long[]} backing is used. */
    public static long MAX_SMALL_CAPACITY = Integer.MAX_VALUE - 8;

    /**
     * Handle used to atomically OR a bit into a single backing word, so that {@link
     * #setAndTestWasUnset(long)} can run lock-free from multiple threads without losing concurrent
     * updates.
     */
    private static final VarHandle LONG_ARRAY_HANDLE = MethodHandles.arrayElementVarHandle(long[].class);

    // All made public for inlining (optimization):
    /** The current capacity in 64-bit words. */
    public long size;
    /** The small backing array, or {@code null} when the large backing is used. */
    public long[] bits;
    /** The large (big-array) backing, or {@code null} when the small backing is used. */
    public long[][] largeBits;

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
        ensureCapacity(initialSize, enforceLarge);
    }

    /**
     * Sets all bits to zero.
     */
    public void clear() {
        if (largeBits != null) {
            BigArrays.fill(largeBits, 0L);
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
            size = newSize;
            if (size > MAX_SMALL_CAPACITY || enforceLarge || largeBits != null) {
                if (bits != null) {
                    largeBits = BigArrays.wrap(bits);
                }
                bits = null;
                largeBits = BigArrays.ensureCapacity(largeBits == null ? LongBigArrays.EMPTY_BIG_ARRAY : largeBits,
                        size);
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
        if (largeBits != null) {
            long arrayIndex = index >>> 6;
            BigArrays.set(largeBits, arrayIndex, BigArrays.get(largeBits, arrayIndex) & ~(1L << (index & 0b111111)));
        } else {
            int arrayIndex = (int) (index >>> 6);
            bits[arrayIndex] = bits[arrayIndex] & ~(1L << (index & 0b111111));
        }
    }

    /**
     * Sets (to 1) the bit at the given index.
     *
     * @param index the bit index to set
     */
    public final void set(final long index) {
        if (largeBits != null) {
            long arrayIndex = index >>> 6;
            BigArrays.set(largeBits, arrayIndex, BigArrays.get(largeBits, arrayIndex) | (1L << (index & 0b111111)));
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
     * Atomically sets (to 1) the bit at the given index and reports whether this call was the one
     * that set it, i.e. whether the bit was previously 0. Unlike {@link #set(long)} this is safe for
     * concurrent use by multiple threads: the word is updated with an atomic OR so no concurrent bit
     * set on the same word is ever lost.
     *
     * @param index the bit index to set
     * @return {@code true} if the bit transitioned from 0 to 1, {@code false} if it was already set
     */
    public final boolean setAndTestWasUnset(final long index) {
        final long bit = 1L << (index & 0b111111);
        if (largeBits != null) {
            final long arrayIndex = index >>> 6;
            final long[] segment = largeBits[BigArrays.segment(arrayIndex)];
            final int displacement = BigArrays.displacement(arrayIndex);
            final long old = (long) LONG_ARRAY_HANDLE.getAndBitwiseOr(segment, displacement, bit);
            return (old & bit) == 0;
        } else {
            final int arrayIndex = (int) (index >>> 6);
            final long old = (long) LONG_ARRAY_HANDLE.getAndBitwiseOr(bits, arrayIndex, bit);
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
        if (largeBits != null) {
            long arrayIndex = index >>> 6;
            return ((BigArrays.get(largeBits, arrayIndex) >> (index & 0b111111)) & 1L) == 1;
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
