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

import java.io.IOException;
import java.io.ObjectOutputStream;
import java.lang.invoke.MethodHandles;
import java.lang.invoke.VarHandle;
import java.util.Arrays;
import java.util.Random;

/**
 * A blocked Bloom filter for k-mers, derived from FastFilter's {@code BlockedBloom} and tuned for
 * lookup speed: every key sets/tests a few bits within a small block of adjacent {@code long} words,
 * so a query touches only one or two cache lines. Uses small ({@code int}-indexed) storage up to
 * {@link #MAX_SMALL_CAPACITY} and fastutil {@link BigArrays} beyond it.
 * <p></p>
 * This implementation is derived from
 * <a href="https://raw.githubusercontent.com/FastFilter/fastfilter_java/refs/heads/master/fastfilter/src/main/java/org/fastfilter/bloom/BlockedBloom.java">BlockedBloom.java</a>
 * The corresponding GitHub project is <a href="https://github.com/FastFilter/fastfilter_java">jastfilter_java</a>
 * is under Apache 2.0 license. It is highly optimized for best classification performance...
 */
public class BlockedKMerBloomFilter extends AbstractCountingKMerProbFilter {
    /** Default false-positive probability. */
    public static double DEFAULT_FPP = 0.01d;
    /** Default number of bits allocated per key. */
    public static int DEFAULT_BITS_PER_KEY = 10;

    private static final long serialVersionUID = 1L;

    /**
     * Handle used to atomically OR bits into a single backing word, so that {@link
     * #putLongIfAbsent(long)} can run lock-free from multiple threads without losing concurrent
     * updates.
     */
    private static final VarHandle LONG_ARRAY_HANDLE = MethodHandles.arrayElementVarHandle(long[].class);

    /** Maximum capacity (in words) that still uses the small {@code int}-indexed storage. */
    public static long MAX_SMALL_CAPACITY = Integer.MAX_VALUE - 8;

    /** Number of bits allocated per key. */
    private final int bitsPerKey;
    /** Hash seed used to derive bit positions. */
    private final long seed;
    /** Number of buckets (words) available for bits. */
    private long buckets;
    /** Small ({@code int}-indexed) bit storage, or {@code null} when large storage is used. */
    private long[] data;
    /** Large ({@code long}-indexed) bit storage, or {@code null} when small storage is used. */
    private long[][] largeData;
    /** Number of keys added via the plain {@link #putLong(long)} path. */
    private long entries;

    /**
     * Creates a filter with {@link #DEFAULT_BITS_PER_KEY} bits per key.
     */
    public BlockedKMerBloomFilter() {
        this(DEFAULT_BITS_PER_KEY);
    }

    /**
     * Creates a filter with the given number of bits per key and a fixed default seed.
     *
     * @param bitsPerKey the number of bits allocated per key
     */
    public BlockedKMerBloomFilter(int bitsPerKey) {
        this(bitsPerKey, new Random(42).nextLong());
    }

    /**
     * Creates a filter with the given number of bits per key and hash seed.
     *
     * @param bitsPerKey the number of bits allocated per key
     * @param seed       the hash seed used to derive bit positions
     */
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

    /**
     * Adds the key to the filter and reports whether it was newly added, setting its bits with atomic
     * ORs. This combines a {@link #containsLong(long)} check with {@link #putLong(long)} in a single
     * pass and is safe to call concurrently from multiple threads: no concurrent bit set is lost, so
     * false negatives cannot occur and the resulting filter state is identical to a plain {@code if
     * (!containsLong(key)) putLong(key)} sequence.
     * <p>
     * The "newly added" flag is exact under single-threaded use. Under concurrent use two threads
     * inserting the same absent key may both observe it as new, so {@link #getEntries()} may
     * marginally over-count; this never affects the filter's membership answers.
     *
     * @param key the k-mer, encoded as a {@code long}, to add
     * @return {@code true} if the key was not already present, {@code false} otherwise
     */
    @Override
    public boolean putLongIfAbsent(long key) {
        long hash = hash(key);
        long start = reduce(hash);
        hash = hash ^ Long.rotateLeft(hash, 32);
        long m1 = (1L << hash) | (1L << (hash >> 6));
        long m2 = (1L << (hash >> 12)) | (1L << (hash >> 18));
        long oldA;
        long oldB;
        if (data != null) {
            int s = (int) start;
            oldA = (long) LONG_ARRAY_HANDLE.getAndBitwiseOr(data, s, m1);
            oldB = (long) LONG_ARRAY_HANDLE.getAndBitwiseOr(data, s + 1 + (int) (hash >>> 60), m2);
        } else {
            oldA = orLarge(start, m1);
            oldB = orLarge(start + 1 + (hash >>> 60), m2);
        }
        // Present iff both mask sets were already fully set before this insert.
        boolean added = ((oldA & m1) != m1) || ((oldB & m2) != m2);
        if (added) {
            countAdded();
        }
        return added;
    }

    /**
     * Atomically ORs the given mask into the large-storage word at the given index and returns the
     * previous word value.
     *
     * @param index the word index into the large storage
     * @param mask  the bit mask to OR in
     * @return the word value before the OR
     */
    private long orLarge(long index, long mask) {
        long[] segment = largeData[BigArrays.segment(index)];
        int displacement = BigArrays.displacement(index);
        return (long) LONG_ARRAY_HANDLE.getAndBitwiseOr(segment, displacement, mask);
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
        resetConcurrentEntries();
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

    /**
     * Computes the hash of the given key.
     *
     * @param x the key to hash
     * @return the (deliberately trivial) hash of the given key.
     */
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
    /**
     * Reduces a hash value to a valid start bucket index.
     *
     * @param v the hash value to reduce
     * @return the start bucket index in {@code [0, buckets)} for the given hash value.
     */
    protected final long reduce(final long v) {
        return Math.abs(v % buckets);
        // Using optimization instead of '%', see:
        // http://lemire.me/blog/2016/06/27/a-fast-alternative-to-the-modulo-reduction/
        // and Line 34 in https://github.com/FastFilter/fastfilter_java/blob/master/fastfilter/src/main/java/org/fastfilter/utils/Hash.java
        // In general, it would NOT work for largeData because of potential long-overflow due to the multiplication.
        // return (((((int) v) & 0xffffffffL) * (buckets & 0xffffffffL)) >>> 32);
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
        return entries + concurrentEntryCount();
    }

    /**
     * Folds any keys counted via the concurrent insert path into {@link #entries} before the default
     * serialization runs, so that a deserialized filter reports the correct entry count even though
     * the concurrent counter is transient.
     *
     * @param out the stream the filter is written to
     * @throws IOException if writing fails
     */
    private void writeObject(ObjectOutputStream out) throws IOException {
        entries += drainConcurrentEntries();
        out.defaultWriteObject();
    }
}
