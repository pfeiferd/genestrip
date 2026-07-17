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

import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.util.Arrays;
import java.util.Random;

/**
 * A blocked Bloom filter for k-mers, derived from FastFilter's {@code BlockedBloom} and tuned for
 * lookup speed: every key sets/tests a few bits within a small block of adjacent {@code long} words,
 * so a query touches only one or two cache lines. Uses small ({@code int}-indexed) storage up to
 * {@link #MAX_SMALL_CAPACITY} and a self-managed array of {@code long[]} buckets beyond it. The
 * buckets serve only to address storage beyond a single {@code long[]}; they are made as few and large
 * as the capacity requires and carry no lock — concurrency is provided by a separate pool of stripe
 * locks (see {@link #putLongIfAbsent(long)}).
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

    /** Maximum capacity (in words) that still uses the small {@code int}-indexed storage. */
    public static long MAX_SMALL_CAPACITY = Integer.MAX_VALUE - 8;

    /** Base-2 logarithm of {@link #BUCKET_SIZE}; a word index is split into bucket and displacement here. */
    private static final int BUCKET_SHIFT = 27;
    /** Number of {@code long} words per bucket of the large backing (all but the last are full). */
    private static final int BUCKET_SIZE = 1 << BUCKET_SHIFT;
    /** Mask selecting the in-bucket displacement of a word index. */
    private static final int BUCKET_MASK = BUCKET_SIZE - 1;

    /** Number of stripe locks used to serialize concurrent word updates. Must be a power of two. */
    private static final int SYNCS = 512;
    /** Mask selecting a stripe lock from a word index (power-of-two count, so a cheap AND replaces a modulo). */
    private static final int SYNCS_MASK = SYNCS - 1;

    /** Number of bits allocated per key. */
    private final int bitsPerKey;
    /** Hash seed used to derive bit positions. */
    private final long seed;
    /** Number of buckets (words) available for bits. */
    private long buckets;
    /** Small ({@code int}-indexed) bit storage, or {@code null} when large storage is used. */
    private long[] data;
    /** Large (bucketed) bit storage, or {@code null} when small storage is used. */
    private long[][] largeData;
    /** Number of keys added via the plain {@link #putLong(long)} path. */
    private long entries;
    /**
     * Stripe locks serializing concurrent {@link #putLongIfAbsent(long)} word updates, keyed by word
     * index. Transient (locks carry no state worth serializing) and rebuilt on construction and after
     * deserialization.
     */
    private transient Object[] syncs;

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
        initSyncs();
    }

    private void initSyncs() {
        syncs = new Object[SYNCS];
        for (int i = 0; i < syncs.length; i++) {
            syncs[i] = new Object();
        }
    }

    /**
     * Rebuilds the transient stripe locks after deserialization.
     *
     * @param in the stream to read from
     * @throws IOException            if reading fails
     * @throws ClassNotFoundException if a serialized class cannot be resolved
     */
    private void readObject(ObjectInputStream in) throws IOException, ClassNotFoundException {
        in.defaultReadObject();
        initSyncs();
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
            // Single-threaded classic path: no locking needed.
            setLargeOr(start, m1);
            setLargeOr(start + 1 + (hash >>> 60), m2);
        }
    }

    /**
     * ORs the given mask into the large-storage word at the given index without locking (single-threaded
     * use only).
     *
     * @param index the word index into the large storage
     * @param mask  the bit mask to OR in
     */
    private void setLargeOr(long index, long mask) {
        largeData[(int) (index >>> BUCKET_SHIFT)][(int) (index & BUCKET_MASK)] |= mask;
    }

    /**
     * Adds the key to the filter and reports whether it was newly added, ORing its bits into the
     * backing words under stripe locks. This combines a {@link #containsLong(long)} check with {@link
     * #putLong(long)} in a single pass and is safe to call concurrently from multiple threads: no
     * concurrent bit set is lost, so false negatives cannot occur and the resulting filter state is
     * identical to a plain {@code if (!containsLong(key)) putLong(key)} sequence.
     * <p>
     * Each of the two affected words is updated under the stripe lock ({@link #syncs}) selected by its
     * own word index, so writes to distinct words run in parallel (up to {@link #SYNCS} at a time)
     * regardless of how few buckets the large storage uses. The two words are locked independently and
     * never held at the same time, so no lock-ordering deadlock is possible.
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
            int s2 = s + 1 + (int) (hash >>> 60);
            oldA = orSmall(s, m1);
            oldB = orSmall(s2, m2);
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
     * ORs the given mask into the small-storage word at the given index under its stripe lock, and
     * returns the previous word value.
     *
     * @param index the word index into the small storage
     * @param mask  the bit mask to OR in
     * @return the word value before the OR
     */
    private long orSmall(int index, long mask) {
        synchronized (syncFor(index)) {
            long old = data[index];
            data[index] = old | mask;
            return old;
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
     * ORs the given mask into the large-storage word at the given index under the stripe lock selected
     * by that word's index, and returns the previous word value. Serializing on the word's stripe keeps
     * concurrent writes to the same word from losing updates while letting writes to other words (and
     * hence, typically, other threads) proceed.
     *
     * @param index the word index into the large storage
     * @param mask  the bit mask to OR in
     * @return the word value before the OR
     */
    private long orLarge(long index, long mask) {
        long[] bucket = largeData[(int) (index >>> BUCKET_SHIFT)];
        int displacement = (int) (index & BUCKET_MASK);
        synchronized (syncFor(index)) {
            long old = bucket[displacement];
            bucket[displacement] = old | mask;
            return old;
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
            a = getLarge(start);
            b = getLarge(start + 1 + (hash >>> 60));
        }
        long m1 = (1L << hash) | (1L << (hash >> 6));
        long m2 = (1L << (hash >> 12)) | (1L << (hash >> 18));
        return ((m1 & a) == m1) && ((m2 & b) == m2);
    }

    /**
     * Reads the large-storage word at the given index.
     *
     * @param index the word index into the large storage
     * @return the word value at that index
     */
    private long getLarge(long index) {
        return largeData[(int) (index >>> BUCKET_SHIFT)][(int) (index & BUCKET_MASK)];
    }

    /**
     * Allocates the large (bucketed) backing for {@code words} words on the fixed {@link #BUCKET_SIZE}
     * grid: every bucket holds {@link #BUCKET_SIZE} words except the last, which holds the remainder.
     *
     * @param words the desired capacity in words
     * @return the freshly allocated (zeroed) large backing
     */
    private static long[][] newLargeGrid(long words) {
        int bucketCount = (int) ((words + BUCKET_SIZE - 1) >>> BUCKET_SHIFT);
        long[][] grid = new long[bucketCount][];
        for (int b = 0; b < bucketCount; b++) {
            long startWord = (long) b << BUCKET_SHIFT;
            grid[b] = new long[(int) Math.min(BUCKET_SIZE, words - startWord)];
        }
        return grid;
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
                largeData = newLargeGrid(buckets + 16 + 1);
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
            for (long[] bucket : largeData) {
                Arrays.fill(bucket, 0L);
            }
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
