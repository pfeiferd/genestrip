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
package org.metagene.genestrip.bloom;

import it.unimi.dsi.fastutil.BigArrays;
import it.unimi.dsi.fastutil.longs.LongBigArrays;

import java.util.Arrays;
import java.util.Random;

// This implementation is derived from
// https://raw.githubusercontent.com/FastFilter/fastfilter_java/refs/heads/master/fastfilter/src/main/java/org/fastfilter/bloom/BlockedBloom.java
// The corresponding project from https://github.com/FastFilter/fastfilter_java
// is under Apache 2.0 license.
//
// It is highly optimized for best classification performance...
public class BlockedKMerBloomFilter implements KMerProbFilter {
    private static final long serialVersionUID = 1L;

    public static long MAX_SMALL_CAPACITY = Integer.MAX_VALUE - 8;

    private int bitsPerKey;
    private final long seed;
    private long buckets;
    private long[] data;
    private long[][] largeData;
    private long entries;

    public BlockedKMerBloomFilter(int bitsPerKey) {
        this(bitsPerKey, new Random(42).nextLong());
    }

    public BlockedKMerBloomFilter(int bitsPerKey, long seed) {
        this.bitsPerKey = bitsPerKey;
        this.seed = seed;
        buckets = 0;
        entries = 0;
    }

    @Override
    public void putLong(long key) {
        entries++;
        long hash = hash(key);
        long start = reduce(hash);
        hash = hash ^ Long.rotateLeft(hash, 32);
        long m1 = (1L << hash) | (1L << (hash >> 6));
        long m2 = (1L << (hash >> 12)) | (1L << (hash >> 18));
        if (data != null) {
            int s = (int) start;
            data[s] |= m1;
            data[s + 1 + (int) (hash >>> 60)] |= m2;
        } else {
            BigArrays.set(largeData, start, BigArrays.get(largeData, start) | m1);
            long s2 = start + 1 + (hash >>> 60);
            BigArrays.set(largeData, s2, BigArrays.get(largeData, s2) | m2);
        }
    }

    @Override
    public boolean containsLong(long key) {
        long hash = seed ^ key; // Super simple inlined hash function.
        long start = reduce(hash);
        hash = hash ^ Long.rotateLeft(hash, 32);
        long a;
        long b;
        if (data != null) {
            int s = (int) start;
            a = data[s];
            b = data[s + 1 + (int) (hash >>> 60)];
        } else {
            a = BigArrays.get(largeData, start);
            b = BigArrays.get(largeData, start + 1 + (hash >>> 60));
        }
        long m1 = (1L << hash) | (1L << (hash >> 6));
        long m2 = (1L << (hash >> 12)) | (1L << (hash >> 18));
        return ((m1 & a) == m1) && ((m2 & b) == m2);
    }

    @Override
    public long ensureExpectedSize(long entryCount, boolean enforceLarge) {
        entryCount = Math.max(1, entryCount);
        long bits = entryCount * bitsPerKey;
        long newSize = (bits + 63) / 64;
        entries = 0;
        if (newSize > buckets) {
            buckets = newSize;
            if (buckets + 16 + 1 > MAX_SMALL_CAPACITY || enforceLarge) {
                data = null;
                largeData = BigArrays.ensureCapacity(largeData == null ? LongBigArrays.EMPTY_BIG_ARRAY : largeData,
                        buckets  + 16 + 1);
            } else {
                largeData = null;
                data = new long[(int) buckets + 16 + 1];
            }
        }
        return bits;
    }

    protected final long hash(long x) {
        /*
        x += seed;
        x = (x ^ (x >>> 33)) * 0xff51afd7ed558ccdL;
        x = (x ^ (x >>> 33)) * 0xc4ceb9fe1a85ec53L;
        x = x ^ (x >>> 33);
        return x;
         */
        return seed ^ x;
    }

    // The super simple hash function from above cannot be combined with Lemire's optimization.
    // It results in a very high fpp...
    // The more complex hash function from Lemire outweighs the cost of the modulo operator used here.
    // That's why we keep it simple and leave it as it is.
    protected final long reduce(final long v) {
        if (largeData == null) {
            return (v < 0 ? -v : v) % buckets;
        } else {
            // Using optimization instead of '%', see:
            // http://lemire.me/blog/2016/06/27/a-fast-alternative-to-the-modulo-reduction/
            // and Line 34 in https://github.com/FastFilter/fastfilter_java/blob/master/fastfilter/src/main/java/org/fastfilter/utils/Hash.java
            // In general, it would NOT work for largeData because of potential long-overflow due to the multiplicaton.
            return (((((int) v) & 0xffffffffL) * (buckets & 0xffffffffL)) >>> 32);
        }
    }

    @Override
    public void clear() {
        if (largeData != null) {
            BigArrays.fill(largeData, 0L);
        } else if (data != null) {
            Arrays.fill(data, 0L);
        }
    }

    @Override
    public long getEntries() {
        return entries;
    }

    @Override
    public double getFpp() {
        return 0.01;
    }
}
