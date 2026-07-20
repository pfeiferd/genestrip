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
				xor ? new XORKMerBloomFilter(entryFpp, sumBucketSizes(bucketSizes))
						: new MurmurKMerBloomFilter(entryFpp, sumBucketSizes(bucketSizes)),
				optimizedFpp);
	}

	/**
	 * Sums the per-bucket capacities, which is the total reserved capacity of the store and hence the
	 * expected-insertion count the pre-filter must be sized for.
	 *
	 * @param bucketSizes the capacity of each radix bucket.
	 * @return the total reserved capacity across all buckets.
	 */
	private static long sumBucketSizes(int[] bucketSizes) {
		long total = 0;
		for (int s : bucketSizes) {
			if (s < 0) {
				throw new IllegalArgumentException("Negative bucket size: " + s);
			}
			total += s;
		}
		return total;
	}

	/**
	 * Creates a radix k-mer store with the given radix width, per-bucket capacities and an explicit
	 * probabilistic pre-filter.
	 *
	 * @param k the k-mer length.
	 * @param radixBits the number of low k-mer bits used as the radix index.
	 * @param bucketSizes the capacity of each radix bucket; length must be {@code 2^radixBits}.
	 * @param initialValues the initial list of values.
	 * @param filter the probabilistic pre-filter to use (already sized for the summed bucket
	 *               capacities), or null for none.
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
			// Fail fast (lock-free) - (probably) already stored. Checking for sure would mean a linear
			// scan. Honouring useFilter lets a bulk loader that has pre-counted distinct k-mers turn the
			// (probabilistic) duplicate check off and fill the reserved buckets losslessly.
			return false;
		}
		if (filter != null) {
			// Atomic combined membership-check-and-insert - thread-safe via the filter's own stripe
			// locks, so no global store lock is needed to keep the filter consistent. It marks the k-mer
			// and reports whether it was new; when the filter is consulted (useFilter) a k-mer that was
			// already present is a (probable) duplicate and is dropped, which also settles the race
			// between the lock-free check above and this insert. When useFilter is off the result is
			// ignored (the filter is only populated), matching the previous behaviour.
			boolean newKmer = filter.putLongIfAbsent(kmer);
			if (useFilter && !newKmer) {
				return false;
			}
		}
		long remaining = remainingOf(kmer);
		// Reserve this k-mer's slot by locking on the radix bucket itself, so inserts into different
		// radix buckets run in parallel - only k-mers sharing a radix (hence the same bucket and its
		// bucketFill counter) contend. This replaces the former global lock on the whole store.
		int pos;
		synchronized (bucket) {
			int fill = bucketFill[radix];
			if (fill >= bucket.length) {
				// Reserved capacity for this radix bucket is exhausted. The per-bucket sizes can be a
				// slight under-estimate (Bloom-filter FPP variance between counting and filling), so
				// the k-mer is dropped instead of failing the whole build, like a full KMerSortedArray.
				return false;
			}
			pos = fill;
			bucketFill[radix] = fill + 1;
		}
		// Resolve the value index with a lock-free read: all values are registered up front (the store's
		// initialValues), so a value not found here is a caller error and fails fast. No global entry
		// counter is kept during the fill: the stored k-mer count is the sum of the per-radix bucketFill
		// counters, materialized into 'entries' by optimize() (see getEntries()).
		int vi = getRegisteredValueIndex(value);
		bucket[pos] = entryOf(vi, remaining);
		return true;
	}

	@Override
	public long getEntries() {
		// 'entries' is not maintained during the concurrent fill (putLong keeps only the per-radix
		// bucketFill counters); until optimize() materializes it, derive the count from those counters.
		return sorted ? entries : sumBucketFill();
	}

	private long sumBucketFill() {
		long sum = 0;
		for (int fill : bucketFill) {
			sum += fill;
		}
		return sum;
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
		// Lock on the bucket itself: entries only move within their bucket, so all accesses that might
		// target the same entry share the bucket lock, while updates to different buckets run in parallel.
		synchronized (bucket) {
			long entry = bucket[pos];
			int vi = (int) (entry >>> remainingBits);
			V oldValue = indexMap[vi];
			V newValue = provider.getUpdateValue(oldValue);
			if (newValue == null) {
				throw new NullPointerException("Null is not allowed as a value.");
			}
			if (newValue != oldValue && !newValue.equals(oldValue)) {
				// Values are all pre-registered before the update phase, so this is a lock-free read.
				int newVi = getRegisteredValueIndex(newValue);
				bucket[pos] = entryOf(newVi, entry & remainingMask);
				return true;
			}
			return false;
		}
	}

	/**
	 * Reusable scratch buffers for {@link #updateBatch(BatchBuffers, KMerStore.UpdateValueProvider)}.
	 * A caller (typically one per reader thread) appends up to {@code capacity} k-mers via
	 * {@link #add(long)} and then flushes the batch; reusing one instance keeps the batched update
	 * allocation-free. A {@code BatchBuffers} instance is not thread-safe and must not be shared
	 * between threads.
	 */
	public static final class BatchBuffers {
		private final long[] kmers;
		private final long[] remaining;
		private final long[][] buckets;
		private final int[] lo;
		private final int[] hi;
		private final int[] pos;
		private int count;

		/**
		 * Creates buffers holding up to {@code capacity} k-mers per batch.
		 *
		 * @param capacity the maximum number of k-mers accumulated before a flush is required
		 */
		public BatchBuffers(int capacity) {
			if (capacity < 1) {
				throw new IllegalArgumentException("capacity must be >= 1, got " + capacity);
			}
			kmers = new long[capacity];
			remaining = new long[capacity];
			buckets = new long[capacity][];
			lo = new int[capacity];
			hi = new int[capacity];
			pos = new int[capacity];
		}

		/**
		 * Appends a k-mer to the batch.
		 *
		 * @param kmer the k-mer to add
		 * @return {@code true} if the batch is now full and must be flushed before adding more
		 */
		public boolean add(long kmer) {
			kmers[count] = kmer;
			return ++count == kmers.length;
		}

		/**
		 * Returns whether the batch currently holds no k-mers.
		 *
		 * @return {@code true} if empty
		 */
		public boolean isEmpty() {
			return count == 0;
		}
	}

	/**
	 * Batched variant of {@link #update(long, KMerStore.UpdateValueProvider)} that processes all
	 * k-mers accumulated in {@code b} together. It resolves the radix buckets, queries the pre-filter
	 * and runs the per-bucket binary searches in separate passes over the whole batch, so many
	 * independent cache-missing loads are in flight at the same time (memory-level parallelism)
	 * instead of one at a time as in the per-k-mer {@link #update}. For a provider that is a pure
	 * function of the stored value, the observable effect is identical to calling {@link #update} on
	 * each buffered k-mer in insertion order. The batch is emptied on return.
	 *
	 * @param b        the buffers holding the batch of k-mers to update
	 * @param provider supplies the new value from the currently stored value
	 * @return the number of k-mers whose stored value was changed
	 */
	public int updateBatch(BatchBuffers b, KMerStore.UpdateValueProvider<V> provider) {
		if (!sorted) {
			throw new IllegalStateException("Update only works when optimized.");
		}
		final int n = b.count;
		final long[] kmers = b.kmers;
		final long[] remaining = b.remaining;
		final long[][] buckets = b.buckets;
		final int[] lo = b.lo;
		final int[] hi = b.hi;
		final int[] pos = b.pos;

		// Pass 1: resolve each k-mer's radix bucket and pre-filter it. The radixIndex[] and bloom-filter
		// loads are independent across the batch, so they overlap in flight. pos[i] is set to -1 for
		// k-mers that are absent (null bucket or filtered out) and to -2 for the ones still to locate.
		int active = 0;
		for (int i = 0; i < n; i++) {
			final long kmer = kmers[i];
			final int radix = (int) (kmer & radixMask);
			final long[] bucket = radixIndex[radix];
			if (bucket == null || (filter != null && useFilter && !filter.containsLong(kmer))) {
				pos[i] = -1;
				buckets[i] = null;
				continue;
			}
			buckets[i] = bucket;
			remaining[i] = remainingOf(kmer);
			lo[i] = 0;
			hi[i] = bucketFill[radix] - 1;
			pos[i] = -2;
			active++;
		}

		// Pass 2: interleaved binary search. Each round performs at most one comparison per still-active
		// k-mer, issuing that many independent bucket loads together. A k-mer settles to its found
		// position (>= 0) or to "not found" (-1); the pass ends once none remain active.
		while (active > 0) {
			for (int i = 0; i < n; i++) {
				if (pos[i] != -2) {
					continue;
				}
				final int l = lo[i];
				final int h = hi[i];
				if (l > h) {
					pos[i] = -1;
					active--;
					continue;
				}
				final int mid = (l + h) >>> 1;
				final long midRem = buckets[i][mid] & remainingMask;
				final long rem = remaining[i];
				if (midRem < rem) {
					lo[i] = mid + 1;
				} else if (midRem > rem) {
					hi[i] = mid - 1;
				} else {
					pos[i] = mid;
					active--;
				}
			}
		}

		// Pass 3: apply the provider and write back moved entries, locking on the bucket exactly as
		// update() does. Each entry is re-read under the lock, so duplicate k-mers within the batch
		// compose correctly.
		int moved = 0;
		for (int i = 0; i < n; i++) {
			final int p = pos[i];
			if (p < 0) {
				continue;
			}
			final long[] bucket = buckets[i];
			synchronized (bucket) {
				final long entry = bucket[p];
				final int vi = (int) (entry >>> remainingBits);
				final V oldValue = indexMap[vi];
				final V newValue = provider.getUpdateValue(oldValue);
				if (newValue == null) {
					throw new NullPointerException("Null is not allowed as a value.");
				}
				if (newValue != oldValue && !newValue.equals(oldValue)) {
					// Values are all pre-registered before the update phase, so this is a lock-free read.
					final int newVi = getRegisteredValueIndex(newValue);
					bucket[p] = entryOf(newVi, entry & remainingMask);
					moved++;
				}
			}
		}
		b.count = 0;
		return moved;
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
		// capacity-based seed leaves gaps, so this restores a dense [0, entries) numbering. The running
		// offset ends at the total fill, which is exactly the stored k-mer count, so it also
		// materializes 'entries' (not maintained during the concurrent fill; see getEntries()).
		long offset = 0;
		for (int r = 0; r < radixIndex.length; r++) {
			bucketOffset[r] = offset;
			offset += bucketFill[r];
		}
		entries = offset;
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
	public void setIndexAtPosition(long pos, int index) {
		// A global storage position is bucketOffset[radix] + localPos (see getLong()/visit()); recover
		// the owning bucket and rewrite its entry, keeping the remaining k-mer bits and swapping in the
		// new value index.
		int radix = radixForPos(pos);
		long[] bucket = radixIndex[radix];
		int localPos = (int) (pos - bucketOffset[radix]);
		bucket[localPos] = entryOf(index, bucket[localPos] & remainingMask);
		// This entry's value changed, so any cached per-value k-mer counts are stale.
		invalidateNKmersPerTaxid();
	}

	/**
	 * Returns the radix bucket that owns the given global storage position. {@code bucketOffset} is
	 * non-decreasing, so the owning bucket is the rightmost one whose offset does not exceed
	 * {@code pos} (empty buckets share their successor's offset and are skipped naturally).
	 *
	 * @param pos the global storage position, as reported by {@link #getLong(long, long[])} or
	 *            {@link #visit(KMerStore.IndexedKMerStoreVisitor)}
	 * @return the radix of the bucket containing {@code pos}
	 */
	private int radixForPos(long pos) {
		int lo = 0;
		int hi = bucketOffset.length - 1;
		int radix = 0;
		while (lo <= hi) {
			int mid = (lo + hi) >>> 1;
			if (bucketOffset[mid] <= pos) {
				radix = mid;
				lo = mid + 1;
			} else {
				hi = mid - 1;
			}
		}
		return radix;
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
