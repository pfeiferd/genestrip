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

import it.unimi.dsi.fastutil.BigArrays;
import it.unimi.dsi.fastutil.longs.LongBigArrays;

public class LargeBitVector implements Serializable {
    private static final long serialVersionUID = 1L;

    public static long MAX_SMALL_CAPACITY = Integer.MAX_VALUE - 8;

    // All made public for inlining (optimization):
    public long size;
    public long[] bits;
    public long[][] largeBits;

    public LargeBitVector(long initialSize) {
        this(initialSize, false);
    }

    public LargeBitVector(long initialSize, boolean enforceLarge) {
        size = -1;
        ensureCapacity(initialSize, enforceLarge);
    }

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

    public boolean isLarge() {
        return largeBits != null;
    }

    public long getBitSize() {
        return size * 64;
    }

    public void clear(long index) {
        if (largeBits != null) {
            long arrayIndex = index >>> 6;
            BigArrays.set(largeBits, arrayIndex, BigArrays.get(largeBits, arrayIndex) & ~(1L << (index & 0b111111)));
        } else {
            int arrayIndex = (int) (index >>> 6);
            bits[arrayIndex] = bits[arrayIndex] & ~(1L << (index & 0b111111));
        }
    }

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
