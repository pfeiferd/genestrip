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
 * Bloom filter whose every key sets all of its bits within one single {@code long} of the backing
 * array, so that a lookup as well as an insert touches exactly <em>one</em> word, where
 * {@link BlockedKMerBloomFilter} touches two.
 * <p>
 * A key's hash selects its word ({@code hash mod k}) and, from separate hash bits, the
 * {@link #getHashBits() n} bit positions it sets within that word's 64 bits. Insertion ORs those bits
 * in, so bits are only ever added and a key stays findable no matter what else lands in its word -
 * the filter therefore keeps the {@link KMerProbFilter} contract of never producing a false negative.
 * A lookup reports the key as present if all of its bits are set.
 * <p>
 * Confining a key to one word costs some accuracy compared to spreading its bits over the whole
 * array: the words are loaded unevenly, and an above-average loaded word raises the false-positive
 * rate for every key in it. That is the price for the single access.
 * <p>
 * Insertion via {@link #putLong(long)} is safe for concurrent use, locking like
 * {@link BlockedKMerBloomFilter}: the small backing on its single {@code long[]}, the bucketed backing
 * on the bucket owning the affected word. As a key occupies one word only, that lock covers the whole
 * read-modify-write, so the "newly added" flag is exact on both backings. Lookups via
 * {@link #containsLong(long)} stay unsynchronized and may miss a concurrent insert until it is
 * published by other means.
 */
public class SingleWordKMerBloomFilter implements KMerProbFilter {
    /** Default number of bits allocated per key. */
    public static final int DEFAULT_BITS_PER_KEY = 10;
    /**
     * Maximum number of bits a key may set. Every bit position consumes 6 bits of the hash, and the
     * positions are taken from one 64-bit hash, which bounds them at 10.
     */
    public static final int MAX_HASH_BITS = 10;

    private static final long serialVersionUID = 1L;

    /** Fixed hash seed used by the constructors that do not take one. */
    private static final long DEFAULT_SEED = new Random(42).nextLong();

    /** Maximum capacity (in words) that still uses the small {@code int}-indexed storage. */
    public static final long MAX_SMALL_CAPACITY = Integer.MAX_VALUE - 8;

    /** Minimum base-2 logarithm of the large-backing bucket width (in words). */
    public static final int MIN_BUCKET_SHIFT = 2;
    /** Maximum base-2 logarithm of the large-backing bucket width; keeps a bucket int-indexable. */
    public static final int MAX_BUCKET_SHIFT = 27;

    /** Number of bits allocated per key, i.e. the memory sizing. */
    private final int bitsPerKey;
    /** Number of bits a key sets within its word. */
    private final int hashBits;
    /** Hash seed used to derive word index and bit positions. */
    private final long seed;
    /** Base-2 logarithm of the large-backing bucket width (words per bucket). */
    private final int bucketShift;
    /** Mask selecting the in-bucket displacement of a word index; derived from {@link #bucketShift}. */
    private transient int bucketMask;
    /** Number of words available for bits. */
    private long words;
    /**
     * Small ({@code int}-indexed) bit storage, or {@code null} when large storage is used. Doubles as
     * the lock guarding its words in {@link #putLong(long)}, hence {@code final}: the reference must be
     * safely published so that concurrent inserters cannot lock different objects.
     */
    private final long[] data;
    /**
     * Large (bucketed) bit storage, or {@code null} when small storage is used. Each bucket doubles as
     * the lock guarding its own words in {@link #putLong(long)}.
     */
    private final long[][] largeData;

    /**
     * Creates a filter sized for {@code expectedInsertions} k-mers with
     * {@link #DEFAULT_BITS_PER_KEY} bits per key and the number of bits per key that suits that sizing.
     *
     * @param expectedInsertions the expected number of k-mers to be inserted
     */
    public SingleWordKMerBloomFilter(long expectedInsertions) {
        this(expectedInsertions, DEFAULT_BITS_PER_KEY);
    }

    /**
     * Creates a filter sized for {@code expectedInsertions} k-mers with the given number of bits per
     * key, setting {@link #optimalHashBits(int)} bits per key.
     *
     * @param expectedInsertions the expected number of k-mers to be inserted
     * @param bitsPerKey         the number of bits allocated per key
     */
    public SingleWordKMerBloomFilter(long expectedInsertions, int bitsPerKey) {
        this(expectedInsertions, bitsPerKey, optimalHashBits(bitsPerKey));
    }

    /**
     * Creates a filter sized for {@code expectedInsertions} k-mers that sets the given number of bits
     * per key, with a fixed default seed.
     *
     * @param expectedInsertions the expected number of k-mers to be inserted
     * @param bitsPerKey         the number of bits allocated per key
     * @param hashBits           the number of bits a key sets within its word
     */
    public SingleWordKMerBloomFilter(long expectedInsertions, int bitsPerKey, int hashBits) {
        this(expectedInsertions, bitsPerKey, hashBits, DEFAULT_SEED, minBucketShift(expectedInsertions, bitsPerKey),
                false);
    }

    /**
     * Creates a filter with the given sizing, seed and large-backing bucket width. The backing (small
     * {@code int}-indexed array vs. bucketed {@code long[][]}) is chosen from the resulting word count
     * relative to {@link #MAX_SMALL_CAPACITY} and is fixed for the lifetime of the filter.
     *
     * @param expectedInsertions the expected number of k-mers to be inserted
     * @param bitsPerKey         the number of bits allocated per key
     * @param hashBits           the number of bits a key sets within its word
     * @param seed               the hash seed
     * @param bucketShift        base-2 logarithm of the large-backing bucket width in words
     * @param forceLarge         whether to use the bucketed backing even when the small one would suffice
     */
    SingleWordKMerBloomFilter(long expectedInsertions, int bitsPerKey, int hashBits, long seed, int bucketShift,
            boolean forceLarge) {
        if (bucketShift < MIN_BUCKET_SHIFT || bucketShift > MAX_BUCKET_SHIFT) {
            throw new IllegalArgumentException(
                    "bucketShift must be in [" + MIN_BUCKET_SHIFT + ", " + MAX_BUCKET_SHIFT + "], got " + bucketShift);
        }
        if (bitsPerKey < 1) {
            throw new IllegalArgumentException("bitsPerKey must be >= 1, got " + bitsPerKey);
        }
        if (hashBits < 1 || hashBits > MAX_HASH_BITS) {
            throw new IllegalArgumentException("hashBits must be in [1, " + MAX_HASH_BITS + "], got " + hashBits);
        }
        this.bitsPerKey = bitsPerKey;
        this.hashBits = hashBits;
        this.seed = seed;
        this.bucketShift = bucketShift;
        this.bucketMask = (1 << bucketShift) - 1;

        words = requiredWords(expectedInsertions, bitsPerKey);
        if (forceLarge || words > MAX_SMALL_CAPACITY) {
            data = null;
            largeData = newLargeGrid(words);
        } else {
            largeData = null;
            data = new long[(int) words];
        }
    }

    /**
     * Returns the number of bits a key should set for the given sizing.
     * <p>
     * Deliberately <em>not</em> the classical Bloom optimum of {@code bitsPerKey * ln 2}: that formula
     * assumes a key's bits are spread over the whole array, whereas here they crowd into one word, so
     * every additional bit fills that word faster. Measured over sizings from 9 to 20 bits per key the
     * false-positive rate bottoms out at 5 to 6 bits and rises again beyond, while the classical
     * formula would ask for up to 10 - at 10 bits per key it recommends 7, which measures 2.21% where 5
     * bits measure 1.80%. Hence half the sizing, capped at 6.
     *
     * @param bitsPerKey the number of bits allocated per key
     * @return the number of bits a key should set
     */
    public static int optimalHashBits(int bitsPerKey) {
        return Math.max(1, Math.min(6, (int) Math.round(bitsPerKey / 2.0)));
    }

    /**
     * Returns the number of words a filter of the given sizing allocates. Unlike
     * {@link BlockedKMerBloomFilter} no padding is needed, as a key never reaches beyond its word.
     *
     * @param expectedInsertions the expected number of k-mers to be inserted
     * @param bitsPerKey         the number of bits allocated per key
     * @return the word count of the backing storage, at least one
     */
    private static long requiredWords(long expectedInsertions, int bitsPerKey) {
        long entryCount = Math.max(1, expectedInsertions);
        return Math.max(1, (entryCount * bitsPerKey + 63) / 64);
    }

    /**
     * Returns the smallest bucket-width exponent that still holds all the words a filter of the given
     * sizing needs, clamped to {@code [}{@link #MIN_BUCKET_SHIFT}{@code , }{@link #MAX_BUCKET_SHIFT}{@code ]}.
     *
     * @param expectedInsertions the expected number of k-mers to be inserted
     * @param bitsPerKey         the number of bits allocated per key
     * @return the smallest sufficient bucket-width exponent
     */
    public static int minBucketShift(long expectedInsertions, int bitsPerKey) {
        long words = requiredWords(expectedInsertions, bitsPerKey);
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
     * Creates a bucket-backed filter, for tests that need the bucketed backing at a size that fits in
     * memory.
     *
     * @param expectedInsertions the expected number of k-mers to be inserted
     * @param bitsPerKey         the number of bits allocated per key
     * @return a filter of the given sizing that uses the bucketed backing
     */
    static SingleWordKMerBloomFilter newLargeBacked(long expectedInsertions, int bitsPerKey) {
        return new SingleWordKMerBloomFilter(expectedInsertions, bitsPerKey, optimalHashBits(bitsPerKey), DEFAULT_SEED,
                minBucketShift(expectedInsertions, bitsPerKey), true);
    }

    /**
     * Returns whether this filter uses the bucketed backing rather than the small one.
     *
     * @return whether this filter uses the bucketed (large) backing
     */
    boolean isLargeBacked() {
        return largeData != null;
    }

    /**
     * Returns the word at the given index, letting a test check which words a key touched.
     *
     * @param index the word index
     * @return the word's content
     */
    long getWord(long index) {
        return data != null ? data[(int) index]
                : largeData[(int) (index >>> bucketShift)][(int) (index & bucketMask)];
    }

    /**
     * Returns the number of bits allocated per key.
     *
     * @return the number of bits allocated per key
     */
    public int getBitsPerKey() {
        return bitsPerKey;
    }

    /**
     * Returns the number of bits a key sets within its word.
     *
     * @return the number of bits a key sets
     */
    public int getHashBits() {
        return hashBits;
    }

    @Override
    public long getBitSize() {
        return words * 64;
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
     * Adds the key to the filter and reports whether it was newly added, ORing all of its bits into its
     * word. This combines a {@link #containsLong(long)} check with the put in a single pass; the
     * resulting filter state is identical to a plain {@code if (!containsLong(key)) putLong(key)}
     * sequence.
     * <p>
     * <strong>Thread-safe:</strong> the whole read-modify-write of the word runs under the lock of its
     * backing, so no concurrent insert is lost and, since a key touches one word only, the returned
     * flag is exact.
     *
     * @param key the k-mer, encoded as a {@code long}, to add
     * @return {@code true} if the key was not already present, {@code false} otherwise
     */
    @Override
    public boolean putLong(long key) {
        long hash = hash(key);
        long mask = mask(hash);
        if (data != null) {
            int index = reduceInt(hash);
            synchronized (data) {
                long old = data[index];
                if ((old & mask) == mask) {
                    return false;
                }
                data[index] = old | mask;
                return true;
            }
        }
        long index = reduce(hash);
        long[] bucket = largeData[(int) (index >>> bucketShift)];
        int displacement = (int) (index & bucketMask);
        synchronized (bucket) {
            long old = bucket[displacement];
            if ((old & mask) == mask) {
                return false;
            }
            bucket[displacement] = old | mask;
            return true;
        }
    }

    @Override
    public boolean containsLong(long key) {
        long hash = hash(key);
        long word;
        if (data != null) {
            word = data[reduceInt(hash)];
        } else {
            long index = reduce(hash);
            word = largeData[(int) (index >>> bucketShift)][(int) (index & bucketMask)];
        }
        long mask = mask(hash);
        return (word & mask) == mask;
    }

    /**
     * Returns the mask of the {@link #hashBits} bits the given hash's key sets within its word. Every
     * bit position takes 6 bits of the hash, which {@code 1L <<} consumes as the shift distance on its
     * own, so the positions are simply successive 6-bit groups. Positions may coincide, which costs a
     * little accuracy but no correctness - as in any Bloom filter whose hashes collide.
     *
     * @param hash the key's hash
     * @return the mask of the bits the key sets
     */
    private long mask(long hash) {
        // Mixed so that the bit positions do not correlate with the word index, which reduce() takes
        // from the same hash - exactly as BlockedKMerBloomFilter does it.
        long mixed = hash ^ Long.rotateLeft(hash, 32);
        // Unrolled and falling through rather than looping with a running shift: every position is an
        // independent shift of 'mixed', so they issue in parallel instead of forming a dependency chain
        // of hashBits steps. (1L << x uses the low 6 bits of x on its own, so no masking is needed.)
        long mask = 0;
        switch (hashBits) {
        case 10:
            mask |= 1L << (mixed >>> 54);
        case 9:
            mask |= 1L << (mixed >>> 48);
        case 8:
            mask |= 1L << (mixed >>> 42);
        case 7:
            mask |= 1L << (mixed >>> 36);
        case 6:
            mask |= 1L << (mixed >>> 30);
        case 5:
            mask |= 1L << (mixed >>> 24);
        case 4:
            mask |= 1L << (mixed >>> 18);
        case 3:
            mask |= 1L << (mixed >>> 12);
        case 2:
            mask |= 1L << (mixed >>> 6);
        default:
            mask |= 1L << mixed;
        }
        return mask;
    }

    /**
     * Allocates the large (bucketed) backing for {@code words} words.
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
     * Computes the hash of the given key, as trivially as {@link BlockedKMerBloomFilter} does.
     *
     * @param x the key to hash
     * @return the (deliberately trivial) hash of the given key
     */
    protected final long hash(long x) {
        return seed ^ x;
    }

    /**
     * Reduces a hash value to a valid word index on the <em>large</em> backing, using Lemire's
     * multiply-shift reduction widened to 64 bits: the upper half of the 128-bit product of the hash
     * (read as unsigned) with {@link #words}, which {@link Math#multiplyHigh(long, long)} yields, lands
     * uniformly in {@code [0, words)}. So no {@code 2^32} limit applies and the 64-bit division of a
     * modulo is avoided. The hash is mixed by a multiplication first because the reduction is driven by its
     * high bits, whereas the trivial {@code seed ^ x} hash carries the k-mer's entropy in the low ones,
     * which the mixing multiplication carries upwards.
     *
     * @param v the hash value to reduce
     * @return the word index in {@code [0, words)} for the given hash value
     */
    protected final long reduce(final long v) {
        // Multiply-shift is driven by the high bits of its input, so the trivial 'seed ^ key' hash is
        // mixed into them first: a multiplication propagates every input bit upwards through the
        // carries. A k-mer carries its entropy in the low bits, which would otherwise hardly reach the
        // index at all.
        long mixedForIndex = v * 0x9E3779B97F4A7C15L;
        // Math.multiplyHigh is signed, so a negative left operand needs the range added back to reach
        // the unsigned product's upper half; 'words' is always positive, so only that side corrects.
        return Math.multiplyHigh(mixedForIndex, words) + ((mixedForIndex >> 63) & words);
    }

    /**
     * Reduces a hash value to a valid word index on the <em>small</em> ({@code int}-indexed) backing,
     * using Lemire's multiply-shift alternative to the modulo: the top 32 bits of
     * {@code (loWord(v) * words)} land uniformly in {@code [0, words)}. Only valid where
     * {@code words < 2^31} keeps the multiply from overflowing, hence the large backing uses
     * {@link #reduce(long)}, which widens the same idea to 64 bits.
     * <p>
     * <strong>Do not carry this reduction over to {@link AbstractKMerBloomFilter}.</strong> Like every
     * multiply-shift it consumes only part of the hash - here its low 32 bits - and so depends on where
     * a key's entropy sits. {@link XORKMerBloomFilter}'s {@code hashFactors[i] ^ x} does not mix at all,
     * so there a reduction of this family costs orders of magnitude of false-positive rate and, through
     * the store's filter-based deduplication, actual k-mers. That is why
     * {@link AbstractKMerBloomFilter#reduce(long)} is and stays a modulo - see there for the numbers.
     *
     * @param v the hash value to reduce
     * @return the word index in {@code [0, words)} for the given hash value
     */
    protected final int reduceInt(final long v) {
        return (int) (((v & 0xffffffffL) * words) >>> 32);
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
