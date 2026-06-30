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
 * {@value #REMAINING_BITS} bits and the value index in their high 19 bits
 * ({@code entry = (valueIndex << }{@value #REMAINING_BITS}{@code ) | remaining}). The remaining bits
 * fit in {@value #REMAINING_BITS} bits as long as {@code 2*k - radixBits <= }{@value #REMAINING_BITS}
 * (validated by the constructor). This is the role that {@code valueIndexes}/{@code largeValueIndexes}
 * play in {@link KMerSortedArray}; here it lives in the spare 16 bits of each k-mer entry, so no
 * separate value-index array (and no large/small variant) is needed — each radix bucket holds at
 * most {@code entries / 2^radixBits} k-mers and therefore fits in a plain {@code int}-indexed
 * {@code long[]} even for very large databases.
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
 */
public class RadixKMerStore<V extends Serializable> extends AbstractKMerStore<V> {
	/**
	 * Minimum (and default) number of low k-mer bits used as the radix index. At least 17 bits are
	 * used so that the remaining k-mer bits ({@code 2*k - radixBits}, at most {@code 2*31 - 17 = 45}
	 * since {@code k <= 31}) fit in the {@value #REMAINING_BITS} low bits of each entry, leaving the
	 * high 19 bits free for the value index.
	 */
	public static final int MIN_RADIX_BITS = 17;
	/** Default number of (low) k-mer bits used as the radix index. */
	public static final int DEFAULT_RADIX_BITS = MIN_RADIX_BITS;
	/** Number of (remaining) k-mer bits stored in the low bits of each entry. */
	public static final int REMAINING_BITS = 45;
	/** Mask selecting the {@value #REMAINING_BITS} stored (remaining) k-mer bits of an entry. */
	public static final long REMAINING_MASK = (1L << REMAINING_BITS) - 1;
	// Since k <= 31 (2*k <= 62) the remaining bits need at most 2*31 - MIN_RADIX_BITS = 45, which 45
	// holds exactly (no spare bit, so MIN_RADIX_BITS and the constructor's fit check are load-bearing),
	// leaving 64 - REMAINING_BITS = 19 bits of each entry for the value index - eight times as many
	// distinct values as KMerSortedArray's short-indexed store. (The value index packs into the
	// entry's high bits via an unsigned shift; entries are never compared to the -1L "erroneous k-mer"
	// sentinel, which lives in the full-k-mer space and stays safe because a reconstructed valid k-mer
	// occupies at most 2*k <= 62 bits, i.e. is always < 2^62.)
	public static final int MAX_VALUES = 8 * KMerSortedArray.MAX_VALUES;

	private static final long serialVersionUID = 1L;

	// Orders entries by their remaining k-mer bits (the value index in the high bits is ignored).
	// The masked values are in [0, 2^45) and hence never negative, so a signed compare is correct.
	private static final LongComparator REMAINING_COMPARATOR =
			(a, b) -> Long.compare(a & REMAINING_MASK, b & REMAINING_MASK);

	// Number of low k-mer bits used as the radix index, and the corresponding mask. Configurable per
	// store; the number of buckets is 2^radixBits == radixIndex.length.
	private final int radixBits;
	private final int radixMask;

	// radixIndex[r] holds all k-mers whose low radixBits bits equal r, or null if there are none.
	private final long[][] radixIndex;
	// Number of entries actually stored in each bucket (<= radixIndex[r].length).
	private final int[] bucketFill;

	// Maps (radix, localPos) to a global storage position. Allocated once (hence final). The
	// constructor seeds it with the prefix sums of the bucket capacities (valid while the store is
	// being filled), and optimize() rebuilds it from the *actual* per-bucket fills so that
	// partially filled buckets still yield a dense [0, entries) numbering. It is rebuilt in
	// optimize() rather than maintained per putLong() because incrementing the fill of bucket r
	// shifts the offset of every later bucket - O(#buckets) per inserted k-mer - whereas the
	// offsets are only ever consumed after optimize() (matching / unique-counting).
	private final long[] bucketOffset;

	public RadixKMerStore(int k, int radixBits, int[] bucketSizes, double entryFpp, double optimizedFpp,
						  List<V> initialValues, boolean xor) {
		this(k, radixBits, bucketSizes, initialValues,
				xor ? new XORKMerBloomFilter(entryFpp) : new MurmurKMerBloomFilter(entryFpp), optimizedFpp);
	}

	protected RadixKMerStore(int k, int radixBits, int[] bucketSizes, List<V> initialValues, KMerProbFilter filter,
							 double optimizedFpp) {
		super(k, MAX_VALUES, initialValues, filter, optimizedFpp);
		if (radixBits < MIN_RADIX_BITS || radixBits > 30) {
			throw new IllegalArgumentException("radixBits must be in [" + MIN_RADIX_BITS + ", 30], got " + radixBits);
		}
		// The remaining (2*k - radixBits) k-mer bits must fit in the REMAINING_BITS low bits of the
		// entry (the high 19 bits hold the value index).
		if (2 * k - radixBits > REMAINING_BITS) {
			throw new IllegalArgumentException("radixBits=" + radixBits + " is too small for k=" + k
					+ "; it must be at least " + (2 * k - REMAINING_BITS) + " so the remaining k-mer bits fit.");
		}
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

	public <W extends Serializable> RadixKMerStore(RadixKMerStore<W> org, KMerStore.ValueConverter<W, V> converter) {
		super(org, converter);
		radixBits = org.radixBits;
		radixMask = org.radixMask;
		radixIndex = org.radixIndex;
		bucketFill = org.bucketFill;
		// The offsets depend only on the (shared) bucket capacities, so they can be shared too.
		bucketOffset = org.bucketOffset;
	}

	/**
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

	private static long entryOf(int valueIndex, long remaining) {
		return (((long) valueIndex) << REMAINING_BITS) | remaining;
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
				final long midRem = bucket[mid] & REMAINING_MASK;
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
				if ((bucket[i] & REMAINING_MASK) == remaining) {
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
		return indexMap[(int) (bucket[pos] >>> REMAINING_BITS)];
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
			long midRem = bucket[mid] & REMAINING_MASK;
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
			int vi = (int) (entry >>> REMAINING_BITS);
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
				bucket[pos] = entryOf(newVi, entry & REMAINING_MASK);
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
		for (int r = 0; r < radixIndex.length; r++) {
			long[] bucket = radixIndex[r];
			if (bucket == null) {
				continue;
			}
			int fill = bucketFill[r];
			if (fill > 1) {
				LongArrays.quickSort(bucket, 0, fill, REMAINING_COMPARATOR);
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
					filter.putLong(((bucket[i] & REMAINING_MASK) << radixBits) | r);
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
				long kmer = ((entry & REMAINING_MASK) << radixBits) | r;
				visitor.nextValue(this, kmer, (int) (entry >>> REMAINING_BITS), base + i);
			}
		}
	}
}
