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
 * Insertion via {@link #putLong(long)} is safe for concurrent use: the small backing locks on its
 * single {@code long[]}, the bucketed backing locks on the bucket owning the affected word — the same
 * technique {@link org.metagene.genestrip.util.LargeBitVector} uses, where the bucket is both the
 * storage and the lock. On the large path the lock striping therefore follows the bucket width: a
 * smaller {@link #bucketShift} yields more, narrower buckets and hence more independent locks (see
 * {@link #putLong(long)} for what concurrency does to its return value). Lookups via
 * {@link #containsLong(long)} stay unsynchronized and may miss a concurrent insert until it is
 * published by other means.
 * <p>
 * This implementation is derived from
 * <a href="https://raw.githubusercontent.com/FastFilter/fastfilter_java/refs/heads/master/fastfilter/src/main/java/org/fastfilter/bloom/BlockedBloom.java">BlockedBloom.java</a>
 * The corresponding GitHub project is <a href="https://github.com/FastFilter/fastfilter_java">jastfilter_java</a>
 * is under Apache 2.0 license. It is highly optimized for best classification performance...
 */
public class BlockedBloomFilter implements ProbFilter {
    /** Default false-positive probability. */
    public static final double DEFAULT_FPP = 0.01d;
    /** Default number of bits allocated per key. */
    public static final int DEFAULT_BITS_PER_KEY = 10;

    // Bumped from 3: hash() changed from the trivial 'seed ^ x' to the MurmurHash3 finalizer and
    // reduce() dropped the mixing multiplication that used to compensate for it, so a key maps to
    // different words and bits than before.
    //
    // The bump is what was missing when that change was made, and the consequence was severe: a
    // database serialized beforehand still deserialized, because the version had not moved, and its
    // store then answered "not present" for about 98% of the k-mers it actually held - the store
    // pre-filters every lookup with this filter (see RadixKMerStore.locateBatch). Lookups do not fail
    // loudly, they just return nothing, so k-mers went missing silently all the way into the FT
    // intersection counts.
    //
    // Hence: whenever hash(), reduce(), reduceInt() or the bit derivation in putLong() changes, this
    // version must be bumped, so that an incompatible filter fails to load instead of answering
    // wrongly. Filters written by an older version have to be regenerated - as do the databases that
    // embed them.
    private static final long serialVersionUID = 4L;

    /** Fixed hash seed used by the constructors that do not take one. */
    static final long DEFAULT_SEED = new Random(42).nextLong();

    /** Maximum capacity (in words) that still uses the small {@code int}-indexed storage. */
    public static final long MAX_SMALL_CAPACITY = Integer.MAX_VALUE - 8;

    /**
     * Minimum base-2 logarithm of the large-backing bucket width (in {@code long} words). The grid only
     * ever allocates as many buckets as the requested capacity needs (see {@link #newLargeGrid(long)}),
     * so a small width simply yields more, smaller buckets without wasting memory. It must keep a
     * bucket wider than the {@link #MAX_WORD_SPAN}-word span a key touches, because both of a key's
     * words are placed in the same bucket (see {@link #secondWord(long[], int, long)}).
     */
    public static final int MIN_BUCKET_SHIFT = 5;

    /**
     * Largest displacement between the two words of a key, i.e. {@code 1 + 15} for the four bits the
     * displacement is taken from. A bucket must hold more words than this so that a key's second word
     * can always be placed in the same bucket as its first.
     */
    private static final int MAX_WORD_SPAN = 16;
    /** Maximum base-2 logarithm of the large-backing bucket width; keeps a bucket int-indexable. */
    public static final int MAX_BUCKET_SHIFT = 27;

    /** Number of bits allocated per key. */
    private final int bitsPerKey;
    /** Hash seed used to derive bit positions. */
    protected final long seed;
    // A large-backing word index splits into bucket (>>> bucketShift) and in-bucket displacement
    // (& bucketMask). bucketShift is fixed at construction; bucketMask is derived (hence transient and
    // rebuilt in readObject).
    /** Base-2 logarithm of the large-backing bucket width (words per bucket). */
    private final int bucketShift;
    /** Mask selecting the in-bucket displacement of a word index; derived from {@link #bucketShift}. */
    private transient int bucketMask;
    /** Number of buckets (words) available for bits. */
    protected long buckets;
    /**
     * Small ({@code int}-indexed) bit storage, or {@code null} when large storage is used. Doubles as
     * the lock guarding its own words in {@link #putLong(long)}, hence {@code final}: the reference must
     * be safely published so that concurrent inserters cannot lock different objects.
     */
    private final long[] data;
    /**
     * Large (bucketed) bit storage, or {@code null} when small storage is used. Each bucket doubles as
     * the lock guarding its own words in {@link #putLong(long)}; the grid is allocated once at
     * construction and never reshaped, so those locks are stable for the filter's lifetime.
     */
    private final long[][] largeData;

    /**
     * Creates a filter sized for {@code expectedInsertions} k-mers with {@link #DEFAULT_BITS_PER_KEY}
     * bits per key.
     *
     * @param expectedInsertions the expected number of k-mers to be inserted
     */
    public BlockedBloomFilter(long expectedInsertions) {
        this(expectedInsertions, DEFAULT_BITS_PER_KEY);
    }

    /**
     * Creates a filter sized for {@code expectedInsertions} k-mers with the given number of bits per
     * key and a fixed default seed.
     *
     * @param expectedInsertions the expected number of k-mers to be inserted
     * @param bitsPerKey the number of bits allocated per key
     */
    public BlockedBloomFilter(long expectedInsertions, int bitsPerKey) {
        this(expectedInsertions, bitsPerKey, DEFAULT_SEED);
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
    public BlockedBloomFilter(long expectedInsertions, int bitsPerKey, long seed) {
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
    public BlockedBloomFilter(long expectedInsertions, int bitsPerKey, long seed, int bucketShift) {
        this(expectedInsertions, bitsPerKey, seed, bucketShift, false);
    }

    /**
     * Creates a filter as {@link #BlockedBloomFilter(long, int, long, int)} does, but able to take
     * the bucketed backing regardless of the filter's size.
     * <p>
     * Reaching that backing through the size alone means allocating more than
     * {@link #MAX_SMALL_CAPACITY} words, i.e. tens of gigabytes, so {@code forceLarge} is what lets a
     * test exercise the bucketed path (and the per-bucket locking that goes with it) at a size that
     * fits in memory. It is deliberately not public: production code selects its backing from the
     * sizing alone.
     *
     * @param expectedInsertions the expected number of k-mers to be inserted
     * @param bitsPerKey  the number of bits allocated per key
     * @param seed        the hash seed used to derive bit positions
     * @param bucketShift base-2 logarithm of the large-backing bucket width in words; must be in
     *                    {@code [}{@link #MIN_BUCKET_SHIFT}{@code , }{@link #MAX_BUCKET_SHIFT}{@code ]}
     * @param forceLarge  whether to use the bucketed backing even when the small one would suffice
     */
    BlockedBloomFilter(long expectedInsertions, int bitsPerKey, long seed, int bucketShift,
                       boolean forceLarge) {
        checkBucketShift(bucketShift);
        this.bitsPerKey = bitsPerKey;
        this.seed = seed;
        this.bucketShift = bucketShift;
        this.bucketMask = (1 << bucketShift) - 1;

        // Clamp to at least one key so the backing (and reduce()'s modulo) never sizes to zero words.
        long entryCount = Math.max(1, expectedInsertions);
        buckets = (entryCount * bitsPerKey + 63) / 64;
        if (forceLarge || buckets + 16 + 1 > MAX_SMALL_CAPACITY) {
            data = null;
            largeData = newLargeGrid(buckets + 16 + 1);
        } else {
            largeData = null;
            data = new long[(int) buckets + 16 + 1];
        }
    }

    /**
     * Creates a bucket-backed filter with the same defaults as
     * {@link #BlockedBloomFilter(long, int)}, for tests that need the bucketed backing at a size
     * that fits in memory.
     *
     * @param expectedInsertions the expected number of k-mers to be inserted
     * @param bitsPerKey the number of bits allocated per key
     * @return a filter of the given sizing that uses the bucketed backing
     */
    static BlockedBloomFilter newLargeBacked(long expectedInsertions, int bitsPerKey) {
        return newLargeBacked(expectedInsertions, bitsPerKey, DEFAULT_SEED,
                minBucketShift(expectedInsertions, bitsPerKey));
    }

    /**
     * Creates a bucket-backed filter with the given seed and bucket width, for tests that need the
     * bucketed backing at a size that fits in memory.
     *
     * @param expectedInsertions the expected number of k-mers to be inserted
     * @param bitsPerKey  the number of bits allocated per key
     * @param seed        the hash seed used to derive bit positions
     * @param bucketShift base-2 logarithm of the large-backing bucket width in words
     * @return a filter of the given sizing that uses the bucketed backing
     */
    static BlockedBloomFilter newLargeBacked(long expectedInsertions, int bitsPerKey, long seed,
                                             int bucketShift) {
        return new BlockedBloomFilter(expectedInsertions, bitsPerKey, seed, bucketShift, true);
    }

    /**
     * Returns whether this filter uses the bucketed backing rather than the small one. Lets a test
     * confirm that it really exercises the bucketed path, which is otherwise indistinguishable from
     * the outside.
     *
     * @return whether this filter uses the bucketed (large) backing
     */
    boolean isLargeBacked() {
        return largeData != null;
    }

    @Override
    public long getBitSize() {
        return buckets * 64;
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

    /**
     * Adds the key to the filter and reports whether it was newly added, ORing its bits into the
     * backing words. This combines a {@link #containsLong(long)} check with put in a
     * single pass; the resulting filter state is identical to a plain
     * {@code if (!containsLong(key)) putLong(key)} sequence.
     * <p>
     * <strong>Thread-safe:</strong> every read-modify-write of a backing word runs under a lock, so no
     * concurrent insert is ever lost, and both of a key's words are covered by one lock: its
     * {@code long[]} on the small backing, the bucket holding both words on the large one. The "newly
     * added" flag is therefore exact on both backings. On the large backing writes to different buckets
     * still proceed in parallel, so the lock striping continues to follow the bucket width.
     *
     * @param key the k-mer, encoded as a {@code long}, to add
     * @return {@code true} if the key was not already present, {@code false} otherwise
     */
    @Override
    public boolean putLong(long key) {
        long hash = hash(key);
        long mixed = hash ^ Long.rotateLeft(hash, 32);
        long m1 = (1L << mixed) | (1L << (mixed >> 6));
        long m2 = (1L << (mixed >> 12)) | (1L << (mixed >> 18));
        long oldA;
        long oldB;
        if (data != null) {
            int s = reduceInt(hash);
            int s2 = s + 1 + (int) (mixed >>> 60);
            synchronized (data) {
                oldA = data[s];
                data[s] = oldA | m1;
                oldB = data[s2];
                data[s2] = oldB | m2;
            }
        } else {
            long start = reduce(hash);
            // One walk through the outer array and one lock for both words, as they share a bucket.
            long[] bucket = largeData[(int) (start >>> bucketShift)];
            int first = (int) (start & bucketMask);
            int second = secondWord(bucket, first, mixed);
            synchronized (bucket) {
                oldA = bucket[first];
                bucket[first] = oldA | m1;
                oldB = bucket[second];
                bucket[second] = oldB | m2;
            }
        }
        // Present iff both mask sets were already fully set before this insert.
        return ((oldA & m1) != m1) || ((oldB & m2) != m2);
    }

    /**
     * Returns the displacement of a key's second word within the bucket holding its first, wrapping
     * around the bucket's end so that both words always share a bucket - which lets
     * {@link #putLong(long)} and {@link #containsLong(long)} reach the outer array once and lock once.
     * The wrap is exact because a bucket holds more than {@link #MAX_WORD_SPAN} words (see
     * {@link #MIN_BUCKET_SHIFT} and {@link #newLargeGrid(long)}), so subtracting its length once always
     * lands back inside it.
     *
     * @param bucket the bucket holding the key's first word
     * @param first  the displacement of the key's first word within that bucket
     * @param mixed  the mixed hash the displacement between the two words is taken from
     * @return the displacement of the key's second word within the same bucket
     */
    private static int secondWord(long[] bucket, int first, long mixed) {
        int second = first + 1 + (int) (mixed >>> 60);
        return second >= bucket.length ? second - bucket.length : second;
    }



    @Override
    public boolean containsLong(long key) {
        long hash = hash(key);
        long mixed = hash ^ Long.rotateLeft(hash, 32);
        long a;
        long b;
        if (data != null) {
            int s = reduceInt(hash);
            a = data[s];
            b = data[s + 1 + (int) (mixed >>> 60)];
        } else {
            long start = reduce(hash);
            // One walk through the outer array for both words, as they share a bucket.
            long[] bucket = largeData[(int) (start >>> bucketShift)];
            int first = (int) (start & bucketMask);
            a = bucket[first];
            b = bucket[secondWord(bucket, first, mixed)];
        }
        long m1 = (1L << mixed) | (1L << (mixed >> 6));
        long m2 = (1L << (mixed >> 12)) | (1L << (mixed >> 18));
        return ((m1 & a) == m1) && ((m2 & b) == m2);
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
            int size = (int) Math.min(bucketSize, words - startWord);
            // The trailing bucket holds only the remainder, which may be shorter than the span between
            // a key's two words; widen it so that secondWord()'s single wrap stays inside it. Costs at
            // most MAX_WORD_SPAN words once, and bucketSize is never smaller than that (MIN_BUCKET_SHIFT).
            grid[b] = new long[Math.max(MAX_WORD_SPAN + 1, size)];
        }
        return grid;
    }

    /**
     * Computes the hash of the given key by seeding it and running it through the MurmurHash3
     * finalizer, which carries a k-mer's low-bit entropy over the whole word. Everything this filter
     * derives from a hash is driven by its high bits - the word index through {@link #reduce(long)} and
     * {@link #reduceInt(long)}, the bit positions through the rotation fold in {@link #putLong(long)} -
     * so the mixing is what makes those sound. {@link XORBlockedBloomFilter} overrides this with a
     * bare exclusive or and consequently reduces by a modulo instead.
     *
     * @param x the key to hash
     * @return the hash of the given key.
     */
    protected long hash(long x) {
        x += seed;
        x = (x ^ (x >>> 33)) * 0xff51afd7ed558ccdL;
        x = (x ^ (x >>> 33)) * 0xc4ceb9fe1a85ec53L;
        x = x ^ (x >>> 33);
        return x;
    }

    /**
     * Reduces a hash value to a valid start bucket index on the <em>large</em> backing, using Lemire's
     * multiply-shift reduction widened to 64 bits: the upper half of the 128-bit product of the hash
     * (read as unsigned) with {@link #buckets} lands uniformly in {@code [0, buckets)}. That upper half
     * is what {@link Math#multiplyHigh(long, long)} yields, so no {@code 2^32} limit applies and the
     * 64-bit division of a modulo is avoided - it dominated this path, which the bucketed backing walks
     * for every lookup and insert.
     * <p>
     * Being a multiply-shift this is driven by the <em>high</em> bits of its input, which is sound only
     * because {@link #hash(long)} mixes a k-mer's low-bit entropy up into them.
     * {@link XORBlockedBloomFilter}, whose hash does not, overrides this with a modulo.
     *
     * @param v the hash value to reduce
     * @return the start bucket index in {@code [0, buckets)} for the given hash value.
     */
    protected long reduce(final long v) {
        // Math.multiplyHigh is signed, so a negative left operand needs the range added back to reach
        // the unsigned product's upper half; 'buckets' is always positive, so only that side corrects.
        return Math.multiplyHigh(v, buckets) + ((v >> 63) & buckets);
    }

    /**
     * Reduces a hash value to a valid start bucket index on the <em>small</em> ({@code int}-indexed)
     * backing, using Lemire's fast alternative to modulo
     * (<a href="http://lemire.me/blog/2016/06/27/a-fast-alternative-to-the-modulo-reduction/">reference</a>):
     * the top 32 bits of {@code (loWord(v) * buckets)} land uniformly in {@code [0, buckets)}. This is
     * only valid for the small backing, where {@code buckets < 2^31}, so the multiply cannot overflow a
     * signed {@code long}; the large path uses {@link #reduce(long)}, which widens the same idea to 64
     * bits.
     * <p>
     * This consumes even less of the hash than {@link #reduce(long)} does - only its low 32 bits - and
     * so depends all the more on {@link #hash(long)} having spread the key's entropy over the whole
     * word. {@link XORBlockedBloomFilter} overrides it with a modulo for that reason.
     *
     * @param v the hash value to reduce
     * @return the start bucket index in {@code [0, buckets)} for the given hash value.
     */
    protected int reduceInt(final long v) {
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
