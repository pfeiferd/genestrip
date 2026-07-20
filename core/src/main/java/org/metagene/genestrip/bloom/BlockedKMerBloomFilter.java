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
import java.util.Arrays;
import java.util.Random;

/**
 * A blocked Bloom filter for k-mers, derived from FastFilter's {@code BlockedBloom} and tuned for
 * lookup speed: every key sets/tests a few bits within a small block of adjacent {@code long} words,
 * so a query touches only one or two cache lines. Uses small ({@code int}-indexed) storage up to
 * {@link #MAX_SMALL_CAPACITY} and a self-managed array of {@code long[]} buckets beyond it. The
 * buckets serve only to address storage beyond a single {@code long[]}; their width defaults to the
 * smallest power of two that holds the whole filter (see {@link #minBucketShift(long, int)}), so the
 * bucketed backing uses as few buckets as the sizing requires, and can also be set explicitly.
 * <p>
 * This filter is not safe for concurrent insertion: unlike {@link AbstractKMerBloomFilter}, its
 * {@link #putLongIfAbsent(long)} does no locking (see that method). It is only ever filled
 * single-threaded (via {@link #putLong(long)} while building an index/database), so no stripe locks
 * are needed.
 * <p></p>
 * This implementation is derived from
 * <a href="https://raw.githubusercontent.com/FastFilter/fastfilter_java/refs/heads/master/fastfilter/src/main/java/org/fastfilter/bloom/BlockedBloom.java">BlockedBloom.java</a>
 * The corresponding GitHub project is <a href="https://github.com/FastFilter/fastfilter_java">jastfilter_java</a>
 * is under Apache 2.0 license. It is highly optimized for best classification performance...
 */
public class BlockedKMerBloomFilter implements KMerProbFilter {
    /** Default false-positive probability. */
    public static double DEFAULT_FPP = 0.01d;
    /** Default number of bits allocated per key. */
    public static int DEFAULT_BITS_PER_KEY = 10;

    private static final long serialVersionUID = 2L;

    /** Maximum capacity (in words) that still uses the small {@code int}-indexed storage. */
    public static long MAX_SMALL_CAPACITY = Integer.MAX_VALUE - 8;

    /**
     * Minimum base-2 logarithm of the large-backing bucket width (in {@code long} words). The grid only
     * ever allocates as many buckets as the requested capacity needs (see {@link #newLargeGrid(long)}),
     * so a small width simply yields more, smaller buckets without wasting memory. It must stay {@code
     * >= 1} because {@code bucketShift == 0} would make {@link #bucketMask} zero, collapsing every word
     * onto one displacement. Widths below the ~17-word span a key touches do cost cache locality on the
     * large path, but that is a performance trade-off, not a correctness one.
     */
    public static final int MIN_BUCKET_SHIFT = 2;
    /** Maximum base-2 logarithm of the large-backing bucket width; keeps a bucket int-indexable. */
    public static final int MAX_BUCKET_SHIFT = 27;

    /** Number of bits allocated per key. */
    private final int bitsPerKey;
    /** Hash seed used to derive bit positions. */
    private final long seed;
    // A large-backing word index splits into bucket (>>> bucketShift) and in-bucket displacement
    // (& bucketMask). bucketShift is fixed at construction; bucketMask is derived (hence transient and
    // rebuilt in readObject).
    /** Base-2 logarithm of the large-backing bucket width (words per bucket). */
    private final int bucketShift;
    /** Mask selecting the in-bucket displacement of a word index; derived from {@link #bucketShift}. */
    private transient int bucketMask;
    /** Number of buckets (words) available for bits. */
    private long buckets;
    /** Small ({@code int}-indexed) bit storage, or {@code null} when large storage is used. */
    private long[] data;
    /** Large (bucketed) bit storage, or {@code null} when small storage is used. */
    private long[][] largeData;

    /**
     * Creates a filter sized for {@code expectedInsertions} k-mers with {@link #DEFAULT_BITS_PER_KEY}
     * bits per key.
     *
     * @param expectedInsertions the expected number of k-mers to be inserted
     */
    public BlockedKMerBloomFilter(long expectedInsertions) {
        this(expectedInsertions, DEFAULT_BITS_PER_KEY);
    }

    /**
     * Creates a filter sized for {@code expectedInsertions} k-mers with the given number of bits per
     * key and a fixed default seed.
     *
     * @param expectedInsertions the expected number of k-mers to be inserted
     * @param bitsPerKey the number of bits allocated per key
     */
    public BlockedKMerBloomFilter(long expectedInsertions, int bitsPerKey) {
        this(expectedInsertions, bitsPerKey, new Random(42).nextLong());
    }

    /**
     * Creates a filter sized for {@code expectedInsertions} k-mers with the given number of bits per
     * key and hash seed, deriving the smallest large-backing bucket width that still holds the whole
     * filter — so the bucketed backing uses as few buckets as the sizing requires (see
     * {@link #minBucketShift(long, int)}).
     *
     * @param expectedInsertions the expected number of k-mers to be inserted
     * @param bitsPerKey the number of bits allocated per key
     * @param seed       the hash seed used to derive bit positions
     */
    public BlockedKMerBloomFilter(long expectedInsertions, int bitsPerKey, long seed) {
        this(expectedInsertions, bitsPerKey, seed, minBucketShift(expectedInsertions, bitsPerKey));
    }

    /**
     * Returns the smallest bucket-width exponent (a power-of-two word count, so the {@code >>>
     * bucketShift} / {@code & bucketMask} addressing keeps working) that still holds all the words a
     * filter of the given sizing needs. This yields the fewest buckets the sizing requires: a single
     * bucket whenever the words fit within one int-addressable block, otherwise the widest permitted
     * bucket ({@link #MAX_BUCKET_SHIFT}) so the grid stays as small as possible. The result is never
     * below {@link #MIN_BUCKET_SHIFT}.
     *
     * @param expectedInsertions the expected number of k-mers to be inserted
     * @param bitsPerKey the number of bits allocated per key
     * @return the smallest sufficient bucket-width exponent, clamped to
     *         {@code [}{@link #MIN_BUCKET_SHIFT}{@code , }{@link #MAX_BUCKET_SHIFT}{@code ]}
     */
    public static int minBucketShift(long expectedInsertions, int bitsPerKey) {
        long words = requiredWords(expectedInsertions, bitsPerKey);
        // ceil(log2(words)) for words >= 2; 0 for words <= 1 (then clamped up to MIN_BUCKET_SHIFT).
        int shift = words <= 1 ? 0 : 64 - Long.numberOfLeadingZeros(words - 1);
        if (shift < MIN_BUCKET_SHIFT) {
            return MIN_BUCKET_SHIFT;
        }
        if (shift > MAX_BUCKET_SHIFT) {
            return MAX_BUCKET_SHIFT;
        }
        return shift;
    }

    /**
     * Returns the number of {@code long} words of bit storage a filter of the given sizing allocates,
     * including the fixed padding a key's bit span may reach into.
     *
     * @param expectedInsertions the expected number of k-mers to be inserted
     * @param bitsPerKey the number of bits allocated per key
     * @return the total word count of the backing storage
     */
    private static long requiredWords(long expectedInsertions, int bitsPerKey) {
        long entryCount = Math.max(1, expectedInsertions);
        return (entryCount * bitsPerKey + 63) / 64 + 16 + 1;
    }

    /**
     * Creates a filter sized for {@code expectedInsertions} k-mers with the given number of bits per
     * key, hash seed and large-backing bucket width. The backing (small {@code int}-indexed array vs.
     * bucketed {@code long[][]}) is chosen here from the resulting word count relative to
     * {@link #MAX_SMALL_CAPACITY}, and is fixed for the lifetime of the filter. The bucket width only
     * affects filters large enough to spill into the bucketed backing; smaller {@code bucketShift}
     * values grow that backing in smaller, more easily allocated steps.
     *
     * @param expectedInsertions the expected number of k-mers to be inserted
     * @param bitsPerKey  the number of bits allocated per key
     * @param seed        the hash seed used to derive bit positions
     * @param bucketShift base-2 logarithm of the large-backing bucket width in words; must be in
     *                    {@code [}{@link #MIN_BUCKET_SHIFT}{@code , }{@link #MAX_BUCKET_SHIFT}{@code ]}
     */
    public BlockedKMerBloomFilter(long expectedInsertions, int bitsPerKey, long seed, int bucketShift) {
        checkBucketShift(bucketShift);
        this.bitsPerKey = bitsPerKey;
        this.seed = seed;
        this.bucketShift = bucketShift;
        this.bucketMask = (1 << bucketShift) - 1;

        // Clamp to at least one key so the backing (and reduce()'s modulo) never sizes to zero words.
        long entryCount = Math.max(1, expectedInsertions);
        buckets = (entryCount * bitsPerKey + 63) / 64;
        if (buckets + 16 + 1 > MAX_SMALL_CAPACITY) {
            data = null;
            largeData = newLargeGrid(buckets + 16 + 1);
        } else {
            largeData = null;
            data = new long[(int) buckets + 16 + 1];
        }
    }

    private static void checkBucketShift(int bucketShift) {
        if (bucketShift < MIN_BUCKET_SHIFT || bucketShift > MAX_BUCKET_SHIFT) {
            throw new IllegalArgumentException(
                    "bucketShift must be in [" + MIN_BUCKET_SHIFT + ", " + MAX_BUCKET_SHIFT + "], got " + bucketShift);
        }
    }

    /**
     * Restores the filter, deriving the transient {@link #bucketMask} from the deserialized
     * {@link #bucketShift}.
     *
     * @param in the stream to read from
     * @throws IOException            if reading fails
     * @throws ClassNotFoundException if a serialized class cannot be resolved
     */
    private void readObject(ObjectInputStream in) throws IOException, ClassNotFoundException {
        in.defaultReadObject();
        bucketMask = (1 << bucketShift) - 1;
    }

    @Override
    public void putLong(long key) {
        long hash = hash(key);
        long mixed = hash ^ Long.rotateLeft(hash, 32);
        long m1 = (1L << mixed) | (1L << (mixed >> 6));
        long m2 = (1L << (mixed >> 12)) | (1L << (mixed >> 18));
        if (data != null) {
            int s = reduceInt(hash);
            data[s] |= m1;
            data[s + 1 + (int) (mixed >>> 60)] |= m2;
        } else {
            // Single-threaded classic path: no locking needed.
            long start = reduce(hash);
            setLargeOr(start, m1);
            setLargeOr(start + 1 + (mixed >>> 60), m2);
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
        largeData[(int) (index >>> bucketShift)][(int) (index & bucketMask)] |= mask;
    }

    /**
     * Adds the key to the filter and reports whether it was newly added, ORing its bits into the
     * backing words. This combines a {@link #containsLong(long)} check with {@link #putLong(long)} in a
     * single pass; the resulting filter state is identical to a plain
     * {@code if (!containsLong(key)) putLong(key)} sequence.
     * <p>
     * <strong>Not thread-safe:</strong> the bit ORs are unsynchronized, so concurrent inserts may lose
     * updates. This filter is only ever filled single-threaded; a filter that needs concurrent
     * inserts uses {@link AbstractKMerBloomFilter} instead (whose bit vector locks per bucket). The
     * "newly added" flag is exact under this single-threaded use.
     *
     * @param key the k-mer, encoded as a {@code long}, to add
     * @return {@code true} if the key was not already present, {@code false} otherwise
     */
    @Override
    public boolean putLongIfAbsent(long key) {
        long hash = hash(key);
        long mixed = hash ^ Long.rotateLeft(hash, 32);
        long m1 = (1L << mixed) | (1L << (mixed >> 6));
        long m2 = (1L << (mixed >> 12)) | (1L << (mixed >> 18));
        long oldA;
        long oldB;
        if (data != null) {
            int s = reduceInt(hash);
            int s2 = s + 1 + (int) (mixed >>> 60);
            oldA = orSmall(s, m1);
            oldB = orSmall(s2, m2);
        } else {
            long start = reduce(hash);
            oldA = orLarge(start, m1);
            oldB = orLarge(start + 1 + (mixed >>> 60), m2);
        }
        // Present iff both mask sets were already fully set before this insert.
        return ((oldA & m1) != m1) || ((oldB & m2) != m2);
    }

    /**
     * ORs the given mask into the small-storage word at the given index and returns the previous word
     * value (single-threaded use only).
     *
     * @param index the word index into the small storage
     * @param mask  the bit mask to OR in
     * @return the word value before the OR
     */
    private long orSmall(int index, long mask) {
        long old = data[index];
        data[index] = old | mask;
        return old;
    }

    /**
     * ORs the given mask into the large-storage word at the given index and returns the previous word
     * value (single-threaded use only).
     *
     * @param index the word index into the large storage
     * @param mask  the bit mask to OR in
     * @return the word value before the OR
     */
    private long orLarge(long index, long mask) {
        long[] bucket = largeData[(int) (index >>> bucketShift)];
        int displacement = (int) (index & bucketMask);
        long old = bucket[displacement];
        bucket[displacement] = old | mask;
        return old;
    }

    @Override
    public boolean containsLong(long key) {
        long hash = seed ^ key; // Super simple inlined hash function.
        long mixed = hash ^ Long.rotateLeft(hash, 32);
        long a;
        long b;
        if (data != null) {
            int s = reduceInt(hash);
            a = data[s];
            b = data[s + 1 + (int) (mixed >>> 60)];
        } else {
            long start = reduce(hash);
            a = getLarge(start);
            b = getLarge(start + 1 + (mixed >>> 60));
        }
        long m1 = (1L << mixed) | (1L << (mixed >> 6));
        long m2 = (1L << (mixed >> 12)) | (1L << (mixed >> 18));
        return ((m1 & a) == m1) && ((m2 & b) == m2);
    }

    /**
     * Reads the large-storage word at the given index.
     *
     * @param index the word index into the large storage
     * @return the word value at that index
     */
    private long getLarge(long index) {
        return largeData[(int) (index >>> bucketShift)][(int) (index & bucketMask)];
    }

    /**
     * Allocates the large (bucketed) backing for {@code words} words on this filter's bucket grid: every
     * bucket holds {@code 1 << bucketShift} words except the last, which holds the remainder.
     *
     * @param words the desired capacity in words
     * @return the freshly allocated (zeroed) large backing
     */
    private long[][] newLargeGrid(long words) {
        int bucketSize = 1 << bucketShift;
        int bucketCount = (int) ((words + bucketSize - 1) >>> bucketShift);
        long[][] grid = new long[bucketCount][];
        for (int b = 0; b < bucketCount; b++) {
            long startWord = (long) b << bucketShift;
            grid[b] = new long[(int) Math.min(bucketSize, words - startWord)];
        }
        return grid;
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

    /**
     * Reduces a hash value to a valid start bucket index on the <em>large</em> backing. Uses a plain
     * modulo because there {@link #buckets} can exceed {@code 2^32}, which would overflow the 32-bit
     * multiply that {@link #reduceInt(long)} relies on. Historically this modulo was also used for the
     * small backing: with the deliberately trivial {@code seed ^ x} hash, Lemire's multiply-shift
     * reduction distributes the low bits worse than the modulo and can raise the false-positive rate,
     * so switching the small path to {@link #reduceInt(long)} trades a little of that margin for speed.
     *
     * @param v the hash value to reduce
     * @return the start bucket index in {@code [0, buckets)} for the given hash value.
     */
    protected final long reduce(final long v) {
        return Math.abs(v % buckets);
    }

    /**
     * Reduces a hash value to a valid start bucket index on the <em>small</em> ({@code int}-indexed)
     * backing, using Lemire's fast alternative to modulo
     * (<a href="http://lemire.me/blog/2016/06/27/a-fast-alternative-to-the-modulo-reduction/">reference</a>):
     * the top 32 bits of {@code (loWord(v) * buckets)} land uniformly in {@code [0, buckets)}. This is
     * only valid for the small backing, where {@code buckets < 2^31}, so the multiply cannot overflow a
     * signed {@code long}; the large path must keep {@link #reduce(long)} (see there).
     *
     * @param v the hash value to reduce
     * @return the start bucket index in {@code [0, buckets)} for the given hash value.
     */
    protected final int reduceInt(final long v) {
        return (int) (((v & 0xffffffffL) * buckets) >>> 32);
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
}
