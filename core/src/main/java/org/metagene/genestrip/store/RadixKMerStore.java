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
package org.metagene.genestrip.store;

import java.io.Serializable;
import java.util.List;

import org.metagene.genestrip.bloom.KMerProbFilter;
import org.metagene.genestrip.bloom.MurmurKMerBloomFilter;
import org.metagene.genestrip.bloom.XORKMerBloomFilter;

import it.unimi.dsi.fastutil.longs.LongArrays;
import it.unimi.dsi.fastutil.longs.LongComparator;

/**
 * A {@link KMerStore} that indexes the lowest {@code radixBits} bits of every k-mer through a fixed
 * radix table and stores the remaining bits of the k-mer together with its value index inside a
 * {@code long}. The radix width {@code radixBits} is configurable per store (constructor argument;
 * see {@link #DEFAULT_RADIX_BITS}).
 * <p>
 * The k-mer is the lowest {@code 2*k} bits of a {@code long} (see
 * {@link org.metagene.genestrip.util.CGAT}). It is split into:
 * <ul>
 * <li>the low {@code radixBits} bits ({@link #radixOf(long, int)}), which directly index the radix
 *     table, and</li>
 * <li>the higher bits ({@code kmer >>> radixBits}), the "remaining" bits, which are stored per
 *     radix bucket.</li>
 * </ul>
 * Using the <em>low</em> bits as the radix means all {@code radixBits} bits are significant for
 * every {@code k} (the low bits of a k-mer are always populated), so the k-mers spread across the
 * full {@code 2^radixBits} buckets independently of {@code k} — unlike a high-bits radix, whose
 * effective width would shrink with smaller {@code k}.
 * <p>
 * Each bucket is a {@code long[]} whose entries pack the remaining k-mer bits in their low
 * {@code remainingBits} bits and the value index in the high {@code 64 - remainingBits} bits
 * ({@code entry = (valueIndex << remainingBits) | remaining}). The width {@code remainingBits} is
 * derived from {@code radixBits} ({@link #remainingBitsForRadix(int)}): it is sized for the
 * worst-case {@code k} (=31), i.e. {@code 2*31 - radixBits}, so the remaining bits of every valid
 * {@code k <= 31} always fit (no per-{@code k} constructor check is needed). Because a larger
 * {@code radixBits} needs fewer remaining bits, it leaves more high bits for the value index, so the
 * number of distinct values the store can hold ({@link #getMaxValues()} /
 * {@link #maxValuesForRadix(int)}) grows with {@code radixBits}. This packing plays the role that
 * {@code valueIndexes}/{@code largeValueIndexes} play in {@link KMerSortedArray}; here it lives in
 * the spare high bits of each k-mer entry, so no separate value-index array (and no large/small
 * variant) is needed — each radix bucket holds at most {@code entries / 2^radixBits} k-mers and
 * therefore fits in a plain {@code int}-indexed {@code long[]} even for very large databases.
 * <p>
 * A lookup first indexes the radix table: a {@code null} bucket means there is no k-mer with that
 * {@code radixBits}-bit suffix, so {@link #getLong} returns {@code null} immediately — even before
 * the (optional) probabilistic pre-filter is queried. Otherwise it binary-searches the bucket
 * (which {@link #optimize()} has sorted by the remaining bits). The hope is that confining the
 * binary search to a single small bucket — rather than searching one global array — keeps the
 * probed positions in cache and reduces the number of probes.
 * <p>
 * The per-bucket capacities are not grown dynamically; they are supplied to the constructor (the
 * caller is expected to have counted the number of distinct k-mers per radix bucket beforehand,
 * using {@link #radixOf(long, int)} with the same {@code radixBits}).
 * <p>
 * Common functionality (value &lt;-&gt; index mapping, counters, the pre-filter and statistics) is
 * inherited from {@link AbstractKMerStore}.
 *
 * @param <V> the value type mapped to each k-mer
 */
public class RadixKMerStore<V extends Serializable> extends AbstractKMerStore<V> {
	/**
	 * Minimum number of low k-mer bits used as the radix index. The lower this bound, the wider the
	 * reserved remaining-bits field (sized for the worst-case {@code k}, see
	 * {@link #remainingBitsForRadix(int)}) and hence the fewer high bits are left for the value
	 * index: at 16 bits the store supports {@code 2^18} distinct values (see
	 * {@link #maxValuesForRadix(int)}).
	 */
	public static final int MIN_RADIX_BITS = 16;
	/** Maximum number of low k-mer bits used as the radix index. */
	public static final int MAX_RADIX_BITS = 30;
	/**
	 * Default number of (low) k-mer bits used as the radix index. Deliberately kept above
	 * {@link #MIN_RADIX_BITS} so the default store retains its larger value capacity
	 * ({@code 2^19} distinct values) rather than the {@code 2^18} of the minimum width.
	 */
	public static final int DEFAULT_RADIX_BITS = 17;

	// Worst-case number of populated k-mer bits (2 * 31, since k <= 31). The remaining-bits field is
	// sized against this so it holds the remaining bits of every valid k regardless of the store's k.
	private static final int MAX_KMER_BITS = 62;
	// Cap on the value-index width so getMaxValues() stays a positive int / valid array length even
	// for the widest radix (a 30-bit index already allows > 10^9 distinct values). Also keeps the
	// 1 << valueBits below the 32-bit shift wrap-around at radixBits == MAX_RADIX_BITS.
	private static final int MAX_VALUE_INDEX_BITS = 30;

	private static final long serialVersionUID = 1L;

	// Number of low k-mer bits used as the radix index, and the corresponding mask. Configurable per
	// store; the number of buckets is 2^radixBits == radixIndex.length.
	/** Number of low k-mer bits used as the radix index. */
	private final int radixBits;
	/** Bit mask selecting the low {@code radixBits} bits used as the radix index. */
	private final int radixMask;
	// Number of low entry bits reserved for the (stored) remaining k-mer bits, and its mask. Derived
	// from radixBits (see remainingBitsForRadix); the value index occupies the entry bits above them.
	// The value index packs into those high bits via an unsigned shift, so an entry may be negative;
	// entries are never compared to the -1L "erroneous k-mer" sentinel, which lives in the full-k-mer
	// space and stays safe because a reconstructed valid k-mer occupies at most 2*k <= 62 bits.
	/** Number of low entry bits reserved for the stored remaining k-mer bits. */
	private final int remainingBits;
	/** Bit mask selecting the low {@code remainingBits} bits of an entry. */
	private final long remainingMask;

	// radixIndex[r] holds all k-mers whose low radixBits bits equal r, or null if there are none.
	/** Per-radix buckets: {@code radixIndex[r]} holds all entries whose radix is {@code r}, or null if none. */
	private final long[][] radixIndex;
	// Number of entries actually stored in each bucket (<= radixIndex[r].length).
	/** Number of entries actually stored in each bucket. */
	private final int[] bucketFill;

	// Maps (radix, localPos) to a global storage position. Allocated once (hence final). The
	// constructor seeds it with the prefix sums of the bucket capacities (valid while the store is
	// being filled), and optimize() rebuilds it from the *actual* per-bucket fills so that
	// partially filled buckets still yield a dense [0, entries) numbering. It is rebuilt in
	// optimize() rather than maintained per putLong() because incrementing the fill of bucket r
	// shifts the offset of every later bucket - O(#buckets) per inserted k-mer - whereas the
	// offsets are only ever consumed after optimize() (matching / unique-counting).
	/** Maps each radix bucket to the global storage offset of its first entry. */
	private final long[] bucketOffset;

	/**
	 * Number of low entry bits reserved for the (stored) remaining k-mer bits of a store with the
	 * given radix width. Sized for the worst-case {@code k} (=31, i.e. {@code 2*k = 62}): the
	 * remaining bits are {@code 2*k - radixBits}, maximal at {@code k = 31}, so reserving
	 * {@code 62 - radixBits} is enough for every {@code k <= 31}.
	 *
	 * @param radixBits the radix width of the store.
	 * @return the reserved remaining-bits width, {@code 62 - radixBits}.
	 */
	public static int remainingBitsForRadix(int radixBits) {
		return MAX_KMER_BITS - radixBits;
	}

	/**
	 * Maximum number of distinct values a store with the given radix width can hold. The value index
	 * occupies the entry bits above the reserved remaining bits ({@link #remainingBitsForRadix(int)}),
	 * so a wider radix leaves more bits for it and raises this cap. The width is capped at
	 * {@value #MAX_VALUE_INDEX_BITS} bits so the result stays a positive {@code int} / valid array
	 * length.
	 *
	 * @param radixBits the radix width of the store.
	 * @return {@code 2^valueBits}, where {@code valueBits = min(30, 64 - remainingBitsForRadix(radixBits))}.
	 */
	public static int maxValuesForRadix(int radixBits) {
		checkRadixBits(radixBits);
		int valueBits = Math.min(MAX_VALUE_INDEX_BITS, Long.SIZE - remainingBitsForRadix(radixBits));
		return 1 << valueBits;
	}

	private static void checkRadixBits(int radixBits) {
		if (radixBits < MIN_RADIX_BITS || radixBits > MAX_RADIX_BITS) {
			throw new IllegalArgumentException(
					"radixBits must be in [" + MIN_RADIX_BITS + ", " + MAX_RADIX_BITS + "], got " + radixBits);
		}
	}

	/**
	 * Creates a radix k-mer store with the given radix width and per-bucket capacities, using a
	 * bloom pre-filter with the given false-positive probabilities.
	 *
	 * @param k the k-mer length.
	 * @param radixBits the number of low k-mer bits used as the radix index.
	 * @param bucketSizes the capacity of each radix bucket; length must be {@code 2^radixBits}.
	 * @param entryFpp the target false-positive probability of the entry (pre-)filter.
	 * @param optimizedFpp the target false-positive probability after optimization.
	 * @param initialValues the initial list of values.
	 * @param xor true to use an XOR-based bloom filter, false to use a Murmur-based one.
	 */
	public RadixKMerStore(int k, int radixBits, int[] bucketSizes, double entryFpp, double optimizedFpp,
						  List<V> initialValues, boolean xor) {
		this(k, radixBits, bucketSizes, initialValues,
				xor ? new XORKMerBloomFilter(entryFpp) : new MurmurKMerBloomFilter(entryFpp), optimizedFpp);
	}

	/**
	 * Creates a radix k-mer store with the given radix width, per-bucket capacities and an explicit
	 * probabilistic pre-filter.
	 *
	 * @param k the k-mer length.
	 * @param radixBits the number of low k-mer bits used as the radix index.
	 * @param bucketSizes the capacity of each radix bucket; length must be {@code 2^radixBits}.
	 * @param initialValues the initial list of values.
	 * @param filter the probabilistic pre-filter to use, or null for none.
	 * @param optimizedFpp the target false-positive probability after optimization.
	 */
	protected RadixKMerStore(int k, int radixBits, int[] bucketSizes, List<V> initialValues, KMerProbFilter filter,
							 double optimizedFpp) {
		// maxValuesForRadix validates radixBits (before super allocates), and its result depends on
		// radixBits: a wider radix reserves fewer remaining bits and so admits more distinct values.
		super(k, maxValuesForRadix(radixBits), initialValues, filter, optimizedFpp);
		// The reserved remaining bits are sized for the worst-case k, so 2*k - radixBits always fits;
		// no per-k fit check is needed here.
		this.remainingBits = remainingBitsForRadix(radixBits);
		this.remainingMask = (1L << remainingBits) - 1;
		int radixSize = 1 << radixBits;
		if (bucketSizes.length != radixSize) {
			throw new IllegalArgumentException("bucketSizes must have length 2^radixBits = " + radixSize
					+ ", got " + bucketSizes.length);
		}
		this.radixBits = radixBits;
		this.radixMask = radixSize - 1;
		this.radixIndex = new long[radixSize][];
		this.bucketFill = new int[radixSize];
		this.bucketOffset = new long[radixSize];
		long total = 0;
		for (int r = 0; r < radixSize; r++) {
			int s = bucketSizes[r];
			if (s < 0) {
				throw new IllegalArgumentException("Negative bucket size at radix " + r);
			}
			bucketOffset[r] = total; // start of bucket r = sum of all earlier bucket capacities
			// A zero-size bucket holds no k-mers: keep its long[] null so getLong()/update() can
			// short-circuit on it (a null bucket means "no k-mer with this radix").
			radixIndex[r] = s > 0 ? new long[s] : null;
			total += s;
		}
		this.size = total;
		if (filter != null) {
			filter.ensureExpectedSize(total, false);
			filter.clear();
		}
	}

	/**
	 * Creates a store that shares the radix index of another store while converting its values.
	 *
	 * @param <W> the source store's value type.
	 * @param org the store to copy the radix structure from.
	 * @param converter the converter from the source value type to this store's value type.
	 */
	public <W extends Serializable> RadixKMerStore(RadixKMerStore<W> org, KMerStore.ValueConverter<W, V> converter) {
		super(org, converter);
		radixBits = org.radixBits;
		radixMask = org.radixMask;
		remainingBits = org.remainingBits;
		remainingMask = org.remainingMask;
		radixIndex = org.radixIndex;
		bucketFill = org.bucketFill;
		// The offsets depend only on the (shared) bucket capacities, so they can be shared too.
		bucketOffset = org.bucketOffset;
	}

	/**
	 * Returns the number of low k-mer bits this store uses as its radix index.
	 *
	 * @return the number of low k-mer bits this store uses as its radix index.
	 */
	public int getRadixBits() {
		return radixBits;
	}

	@Override
	public <W extends Serializable> KMerStore<W> convertValues(KMerStore.ValueConverter<V, W> converter) {
		return new RadixKMerStore<>(this, converter);
	}

	@Override
	public void initSize(long size) {
		// The capacity is determined by the per-bucket sizes passed to the constructor; this only
		// (re)prepares the optional pre-filter so that the store behaves like other KMerStores.
		if (size > this.size) {
			throw new IllegalArgumentException("Requested size " + size + " exceeds reserved radix capacity " + this.size + ".");
		}
		if (filter != null) {
			filter.ensureExpectedSize(this.size, false);
			filter.clear();
		}
	}

	// --- Insertion / lookup ---------------------------------------------------

	/**
	 * Returns the radix bucket index of the given k-mer for a store using the given radix width.
	 *
	 * @param kmer the k-mer whose radix bucket index is computed.
	 * @param radixBits the radix width of the store.
	 * @return the radix bucket index of the given k-mer for a store using {@code radixBits} radix
	 *         bits (its low {@code radixBits} bits), in {@code [0, 2^radixBits)}. Callers that
	 *         pre-count k-mers per bucket (to size this store) must use this exact mapping with the
	 *         same {@code radixBits} the store will be created with.
	 */
	public static int radixOf(long kmer, int radixBits) {
		return (int) (kmer & ((1 << radixBits) - 1));
	}

	// The k-mer bits above the radix, i.e. what is stored (sorted) inside a bucket.
	private long remainingOf(long kmer) {
		return kmer >>> radixBits;
	}

	private long entryOf(int valueIndex, long remaining) {
		return (((long) valueIndex) << remainingBits) | remaining;
	}

	@Override
	public boolean putLong(long kmer, V value) {
		if (value == null) {
			throw new NullPointerException("null is not allowed as a value.");
		}
		sorted = false;
		int radix = (int) (kmer & radixMask);
		long[] bucket = radixIndex[radix];
		if (bucket == null) {
			// No capacity was reserved for this radix prefix. With exact per-bucket sizing this
			// never happens; when the sizes are an estimate (e.g. counted via a Bloom filter) a
			// k-mer may still arrive here, in which case it is dropped rather than failing.
			return false;
		}
		if (filter != null && useFilter && filter.containsLong(kmer)) {
			// Fail fast - (probably) already stored. Checking for sure would mean a linear scan.
			// Honouring useFilter lets a bulk loader that has pre-counted distinct k-mers turn the
			// (probabilistic) duplicate check off and fill the reserved buckets losslessly.
			return false;
		}
		long remaining = remainingOf(kmer);
		int pos;
		int vi;
		synchronized (this) {
			if (filter != null && useFilter && filter.containsLong(kmer)) {
				return false;
			}
			int fill = bucketFill[radix];
			if (fill >= bucket.length) {
				// Reserved capacity for this radix bucket is exhausted. The per-bucket sizes can be a
				// slight under-estimate (Bloom-filter FPP variance between counting and filling), so
				// the k-mer is dropped instead of failing the whole build, like a full KMerSortedArray.
				return false;
			}
			pos = fill;
			bucketFill[radix] = fill + 1;
			vi = getAddValueIndex(value);
			if (filter != null) {
				filter.putLong(kmer);
			}
			entries++;
		}
		bucket[pos] = entryOf(vi, remaining);
		return true;
	}

	// Made final for potential (automated) inlining by JVM
	@Override
	public final V getLong(final long kmer, final long[] posStore) {
		final int radix = (int) (kmer & radixMask);
		final long[] bucket = radixIndex[radix];
		if (bucket == null) {
			// No k-mer with this radix prefix - return early, even before the filter.
			return null;
		}
		if (filter != null && useFilter && !filter.containsLong(kmer)) {
			return null;
		}
		final long remaining = remainingOf(kmer);
		final int fill = bucketFill[radix];
		int pos = -1;
		if (sorted) {
			int lo = 0;
			int hi = fill - 1;
			while (lo <= hi) {
				final int mid = (lo + hi) >>> 1;
				final long midRem = bucket[mid] & remainingMask;
				if (midRem < remaining) {
					lo = mid + 1;
				} else if (midRem > remaining) {
					hi = mid - 1;
				} else {
					pos = mid;
					break;
				}
			}
		} else {
			for (int i = 0; i < fill; i++) {
				if ((bucket[i] & remainingMask) == remaining) {
					pos = i;
					break;
				}
			}
		}
		if (pos < 0) {
			return null;
		}
		if (posStore != null) {
			posStore[0] = bucketOffset[radix] + pos;
		}
		return indexMap[(int) (bucket[pos] >>> remainingBits)];
	}

	@Override
	public boolean update(long kmer, KMerStore.UpdateValueProvider<V> provider) {
		if (!sorted) {
			throw new IllegalStateException("Update only works when optimized.");
		}
		int radix = (int) (kmer & radixMask);
		long[] bucket = radixIndex[radix];
		if (bucket == null) {
			return false;
		}
		if (filter != null && useFilter && !filter.containsLong(kmer)) {
			return false;
		}
		long remaining = remainingOf(kmer);
		int fill = bucketFill[radix];
		int pos = -1;
		int lo = 0;
		int hi = fill - 1;
		while (lo <= hi) {
			int mid = (lo + hi) >>> 1;
			long midRem = bucket[mid] & remainingMask;
			if (midRem < remaining) {
				lo = mid + 1;
			} else if (midRem > remaining) {
				hi = mid - 1;
			} else {
				pos = mid;
				break;
			}
		}
		if (pos < 0) {
			return false;
		}
		// Synchronize only on accesses that (might) target the same entry; see KMerSortedArray.
		synchronized (syncs[((radix * 31) + pos) & (syncs.length - 1)]) {
			long entry = bucket[pos];
			int vi = (int) (entry >>> remainingBits);
			V oldValue = indexMap[vi];
			V newValue = provider.getUpdateValue(oldValue);
			if (newValue == null) {
				throw new NullPointerException("Null is not allowed as a value.");
			}
			if (newValue != oldValue && !newValue.equals(oldValue)) {
				int newVi;
				synchronized (valueMap) {
					newVi = getAddValueIndex(newValue);
					kmersMoved++;
				}
				bucket[pos] = entryOf(newVi, entry & remainingMask);
				return true;
			}
			return false;
		}
	}

	@Override
	public void optimize() {
		if (sorted) {
			return;
		}
		// Orders entries by their remaining k-mer bits (the value index in the high bits is ignored);
		// the masked values are non-negative, so a signed compare is correct. Built locally (rather
		// than kept as a field) so the store stays serializable - a lambda field would not be.
		final long mask = remainingMask;
		final LongComparator remainingComparator = (a, b) -> Long.compare(a & mask, b & mask);
		for (int r = 0; r < radixIndex.length; r++) {
			long[] bucket = radixIndex[r];
			if (bucket == null) {
				continue;
			}
			int fill = bucketFill[r];
			if (fill > 1) {
				LongArrays.quickSort(bucket, 0, fill, remainingComparator);
			}
		}
		sorted = true;
		// Rebuild the position offsets from the actual fills: with partially filled buckets the
		// capacity-based seed leaves gaps, so this restores a dense [0, entries) numbering.
		long offset = 0;
		for (int r = 0; r < radixIndex.length; r++) {
			bucketOffset[r] = offset;
			offset += bucketFill[r];
		}
		// Rework the bloom filter from the fill-time one into the optimized one.
		filter = createOptimizedFilter();
		if (filter != null) {
			for (int r = 0; r < radixIndex.length; r++) {
				long[] bucket = radixIndex[r];
				if (bucket == null) {
					continue;
				}
				int fill = bucketFill[r];
				for (int i = 0; i < fill; i++) {
					// Reassemble the full k-mer: remaining bits shifted back above the radix bits.
					filter.putLong(((bucket[i] & remainingMask) << radixBits) | r);
				}
			}
		}
	}

	@Override
	public void visit(KMerStore.IndexedKMerStoreVisitor<V> visitor) {
		for (int r = 0; r < radixIndex.length; r++) {
			long[] bucket = radixIndex[r];
			if (bucket == null) {
				continue;
			}
			int fill = bucketFill[r];
			long base = bucketOffset[r];
			for (int i = 0; i < fill; i++) {
				long entry = bucket[i];
				// Reassemble the full k-mer: remaining bits shifted back above the radix bits.
				long kmer = ((entry & remainingMask) << radixBits) | r;
				visitor.nextValue(this, kmer, (int) (entry >>> remainingBits), base + i);
			}
		}
	}
}
