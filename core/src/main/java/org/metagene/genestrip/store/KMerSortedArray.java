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

import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.Serializable;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;

import org.metagene.genestrip.bloom.BlockedKMerBloomFilter;
import org.metagene.genestrip.bloom.KMerProbFilter;
import org.metagene.genestrip.bloom.MurmurKMerBloomFilter;
import org.metagene.genestrip.bloom.XORKMerBloomFilter;
import org.metagene.genestrip.util.CGAT;
import org.metagene.genestrip.util.CGATLongBuffer;

import it.unimi.dsi.fastutil.BigArrays;
import it.unimi.dsi.fastutil.BigSwapper;
import it.unimi.dsi.fastutil.longs.LongBigArrays;
import it.unimi.dsi.fastutil.longs.LongComparator;
import it.unimi.dsi.fastutil.objects.Object2LongMap;
import it.unimi.dsi.fastutil.objects.Object2LongOpenHashMap;
import it.unimi.dsi.fastutil.objects.Object2ShortMap;
import it.unimi.dsi.fastutil.objects.Object2ShortOpenHashMap;

public class KMerSortedArray<V extends Serializable> implements KMerStore<V> {
	public static final int MAX_VALUES = ((int) Short.MAX_VALUE) - Short.MIN_VALUE;

	private static final long serialVersionUID = 1L;

	public static long MAX_SMALL_CAPACITY = Integer.MAX_VALUE - 8;

	private long[] kmers;
	private short[] valueIndexes;

	private long[][] largeKmers;
	private short[][] largeValueIndexes;

	private final int k;
	private final Object2ShortMap<V> valueMap;
	// For efficiency and synchronization problems this was turned into an array of
	// lenth Short.MAX_VALUE. Only 0.5 MB in size...
	private final V[] indexMap;
	private final boolean enforceLarge;
	private final double optimizedFpp;

	// Just for optimizing synchronization during updates
	private transient Object[] syncs;

	protected long size;

	private long entries;
	private boolean sorted;
	private boolean initSize;
	private short nextValueIndex;
	private transient KMerProbFilter filter;
	private boolean useFilter;
	// When true, lookups use the radix-guided search instead of plain binary search.
	private boolean radixSearch;
	private Object2LongMap<V> kmerPersTaxid;

	private transient long kmersMoved;

	public KMerSortedArray(int k, double entryFpp, double optimizedFpp, List<V> initialValues, boolean enforceLarge, boolean xor) {
		this(k, initialValues, enforceLarge, xor ? new XORKMerBloomFilter(entryFpp) : new MurmurKMerBloomFilter(entryFpp), optimizedFpp);
	}

	@SuppressWarnings("unchecked")
	protected KMerSortedArray(int k, List<V> initialValues, boolean enforceLarge,
							  KMerProbFilter filter, double optimizedFpp) {
		this.k = k;
		int s = initialValues == null ? 0 : initialValues.size();
		indexMap = (V[]) new Serializable[MAX_VALUES];
		valueMap = new Object2ShortOpenHashMap<>(s);
		initSyncs();
		nextValueIndex = Short.MIN_VALUE;
		this.enforceLarge = enforceLarge;
		if (initialValues != null) {
			for (V v : initialValues) {
				getAddValueIndex(v);
			}
		}
		this.optimizedFpp = optimizedFpp;
		this.filter = filter;
		this.useFilter = filter != null;
		this.kmerPersTaxid = null;
		this.kmersMoved = 0;
	}

	private void readObject(ObjectInputStream ois) throws ClassNotFoundException, IOException {
		ois.defaultReadObject();
		initSyncs();
	}

	private void initSyncs() {
		syncs = new Object[512];
		for (int i = 0; i < syncs.length; i++) {
			syncs[i] = new Object();
		}
	}

	public long getKMersMoved() {
		return kmersMoved;
	}

	@SuppressWarnings("unchecked")
	public <W extends Serializable> KMerSortedArray(KMerSortedArray<W> org, ValueConverter<W, V> converter) {
		kmers = org.kmers;
		valueIndexes = org.valueIndexes;

		largeKmers = org.largeKmers;
		largeValueIndexes = org.largeValueIndexes;

		k = org.k;
		enforceLarge = org.enforceLarge;
		optimizedFpp = org.optimizedFpp;

		size = org.size;

		entries = org.entries;
		sorted = org.sorted;
		initSize = org.initSize;
		nextValueIndex = org.nextValueIndex;
		filter = org.filter;
		useFilter = org.useFilter;
		radixSearch = org.radixSearch;

		indexMap = (V[]) new Serializable[MAX_VALUES];
		valueMap = new Object2ShortOpenHashMap<>(org.valueMap.size());
		syncs = new Object[512];
		for (int i = 0; i < syncs.length; i++) {
			syncs[i] = new Object();
		}

		if (org.kmerPersTaxid != null) {
			kmerPersTaxid = new Object2LongOpenHashMap<>();
		}

		for (short s = 0; s < org.valueMap.size(); s++) {
			W orgValue = org.indexMap[s];
			V value = converter.convertValue(orgValue);
			indexMap[s] = value;
			valueMap.put(value, s);
			if (kmerPersTaxid != null) {
				kmerPersTaxid.put(value, org.kmerPersTaxid.getLong(orgValue));
			}
		}
	}

	public Iterator<V> getValues() {
		return new Iterator<V>() {
			private int i = Short.MIN_VALUE;

			@Override
			public boolean hasNext() {
				return i < nextValueIndex;
			}

			@Override
			public V next() {
				return indexMap[i++ - Short.MIN_VALUE];
			}
		};
	}

	public KMerProbFilter getFilter() {
		return filter;
	}

	public void setFilter(KMerProbFilter filter) {
		this.filter = filter;
	}

	public void setUseFilter(boolean useFilter) {
		this.useFilter = useFilter;
	}

	public boolean isUseFilter() {
		return useFilter;
	}

	/**
	 * Selects the search strategy used by {@link #getLong(long, long[])},
	 * {@link #getLongInlined(long, long[])} and {@link #update(long, UpdateValueProvider)}
	 * for the sorted store.
	 *
	 * @param radixSearch {@code true} for the radix-guided search, {@code false}
	 *                    for the default binary search.
	 */
	public void setRadixSearch(boolean radixSearch) {
		this.radixSearch = radixSearch;
	}

	public boolean isRadixSearch() {
		return radixSearch;
	}

	public Object2LongMap<V> getNKmersPerTaxid() {
		long[] countArray = new long[getNValues()];
		if (largeKmers != null) {
			for (long i = 0; i < entries; i++) {
				countArray[BigArrays.get(largeValueIndexes, i)]++;
			}
		} else {
			for (int i = 0; i < entries; i++) {
				countArray[valueIndexes[i] - Short.MIN_VALUE]++;
			}
		}

		Object2LongMap<V> map = new Object2LongOpenHashMap<>();
		for (int i = 0; i < countArray.length; i++) {
			V value = indexMap[i];
			if (value != null) {
				map.put(value, countArray[i]);
			}
		}
		map.put(null, entries);

		return map;
	}

	public int getNValues() {
		return nextValueIndex - Short.MIN_VALUE;
	}

	public V getValueForIndex(int index) {
		return indexMap[index];
	}

	public int getIndexForValue(V value) {
		short s = valueMap.getOrDefault(value, Short.MAX_VALUE);
		return s == Short.MAX_VALUE ? -1 : (s - Short.MIN_VALUE);
	}

	@Override
	public void initSize(long size) {
		if (initSize) {
			throw new IllegalStateException("Cant initlialize size twice.");
		}
		if (size < 0) {
			throw new IllegalArgumentException("Expected insertions must be > 0.");
		}
		initSize = true;
		filter.ensureExpectedSize(size, false);
		filter.clear();
		this.size = size;

		if (size > MAX_SMALL_CAPACITY || enforceLarge) {
			largeKmers = BigArrays.ensureCapacity(BigArrays.wrap(new long[0]), size);
			largeValueIndexes = BigArrays.ensureCapacity(BigArrays.wrap(new short[0]), size);
		} else {
			kmers = new long[(int) size];
			valueIndexes = new short[(int) size];
		}
	}

	@Override
	public int getK() {
		return k;
	}

	@Override
	public long getEntries() {
		return entries;
	}

	@Override
	public long getSize() {
		return size;
	}

	public boolean isLarge() {
		return largeKmers != null;
	}

	@Override
	public boolean put(CGATLongBuffer buffer, V value) {
		throw new UnsupportedOperationException("Deprecated and so not implemented.");
	}


	@Override
	public boolean put(byte[] nseq, int start, V value) {
		return putLong(CGAT.kMerToLong(nseq, start, k, null), value);
	}

	public boolean isFull() {
		return entries == size;
	}

	/**
	 * @return true if put succeeded, false if (probably) some value is already
	 *         stored under that kmer.
	 */
	public boolean putLong(final long kmer, final V value) {
		if (value == null) {
			throw new NullPointerException("null is not allowed as a value.");
		}
		sorted = false;
		long pos;
		int sindex;
		if (filter.containsLong(kmer)) {
			// Fail fast - we could check if the kmer is indeed stored, but it's way too
			// slow because of linear search in the kmer array...
			return false;
		}
		synchronized (this) {
			// We must check here again for thread safety.
			// It is still better to do it here (again) than not doing it outside to improve
			// parallelism.
			if (filter.containsLong(kmer)) {
				return false;
			}
			if (entries == size) {
				throw new IllegalStateException("Capacity exceeded.");
			}
			pos = entries++;
			sindex = getAddValueIndex(value);
			filter.putLong(kmer);
		}
		if (largeKmers != null) {
			BigArrays.set(largeKmers, pos, kmer);
		} else {
			kmers[(int) pos] = kmer;
		}
		setIndexAtPosition(pos, sindex);
		return true;
	}

	public int getAddValueIndex(V value) {
		short sindex = valueMap.getOrDefault(value, Short.MAX_VALUE);
		if (sindex == Short.MAX_VALUE) {
			if (nextValueIndex == Short.MAX_VALUE) {
				throw new IllegalStateException("Too many different values - only " + MAX_VALUES + " are possible.");
			}
			valueMap.put(value, nextValueIndex);
			indexMap[nextValueIndex - Short.MIN_VALUE] = value;
			sindex = nextValueIndex;
			nextValueIndex++;
		}
		return sindex - Short.MIN_VALUE;
	}

	@Override
	public V get(CGATLongBuffer buffer) {
		throw new UnsupportedOperationException("Deprecated and almost removed.");
	}

	@Override
	public V get(byte[] nseq, int start) {
		return getLong(CGAT.kMerToLong(nseq, start, k, null), null);
	}

	public V get(byte[] nseq, int start, long[] indexStore) {
		return getLong(CGAT.kMerToLong(nseq, start, k, null), indexStore);
	}

	public boolean update(long kmer, UpdateValueProvider<V> provider) {
		if (!sorted) {
			throw new IllegalStateException("Update only works when optimized.");
		}
		if (filter != null && useFilter && !filter.containsLong(kmer)) {
			return false;
		}

		long pos;
		if (largeKmers != null) {
			pos = radixSearch ? radixSearchLarge(kmer, entries)
					: LongBigArrays.binarySearch(largeKmers, 0, entries, kmer);
			if (pos < 0) {
				return false;
			}
		} else {
			pos = radixSearch ? radixSearchSmall(kmer, (int) entries)
					: Arrays.binarySearch(kmers, 0, (int) entries, kmer);
			if (pos < 0) {
				return false;
			}
		}
		// We only must synchronize, if two threads access the same kmer in the same
		// position. But we cannot synchronize on 'pos' because it is a primitive.
		// The helper method getSynchronizationObject() likely returns different objects
		// for different pos values, but always the same object for the same pos value
		// in a multi-threading scenario. This trick greatly decreases synchronization
		// bottlenecks.
		synchronized (syncs[(int) (pos % 512)]) {
			int index;
			if (largeKmers != null) {
				index = BigArrays.get(largeValueIndexes, pos) - Short.MIN_VALUE;
			} else {
				index = valueIndexes[(int) pos] - Short.MIN_VALUE;
			}
			V oldValue = indexMap[index];
			V newValue = provider.getUpdateValue(oldValue);
			if (newValue == null) {
				throw new NullPointerException("Null is not allowed as a value.");
			}
			if (newValue != oldValue && !newValue.equals(oldValue)) {
				// This is important in the multi-threading context:
				synchronized (valueMap) {
					index = getAddValueIndex(newValue);
					kmersMoved++;
				}
				setIndexAtPosition(pos, index);
				return true;
			}
			return false;
		}
	}

	public void setIndexAtPosition(long pos, int index) {
		if (largeKmers != null) {
			BigArrays.set(largeValueIndexes, pos, (short) (index + Short.MIN_VALUE));
		} else {
			valueIndexes[(int) pos] = (short) (index + Short.MIN_VALUE);
		}
	}

	public void fix() {
		fixNKmersPerTaxid();
	}

	protected void fixNKmersPerTaxid() {
		this.kmerPersTaxid = getNKmersPerTaxid();
	}

	public Object2LongMap<V> getFixedNKmersPerTaxid() {
		if (this.kmerPersTaxid == null) {
			fixNKmersPerTaxid();
		}
		return this.kmerPersTaxid;
	}

	public long getKMerAt(long pos) {
		if (largeKmers != null) {
			return BigArrays.get(largeKmers, pos);
		}
		return kmers[(int) pos];
	}

	// Below this interval width the 16-way split degenerates, so the search finishes
	// with plain binary steps.
	private static final int RADIX_MIN_SPAN = 16;

	/**
	 * Radix-guided search over the small (single-array) sorted k-mer store.
	 * <p>
	 * A k-mer occupies the lowest {@code 2*k} bits of the {@code long} (the first base
	 * in the most significant bits, see {@link CGAT#kMerToLongStraight}). Starting from
	 * the top nibble, each step takes the next 4 key bits as a digit {@code d in [0,15]}
	 * and treats the {@code d}-th of 16 equal sub-ranges of the current interval
	 * {@code [lo, hi]} as the predicted location of the key. It probes the start of that
	 * sixteenth ({@code lo + (hi - lo) * d / 16}, the division done with a {@code >>> 4}
	 * shift) and the start of the next one. The two 3-way comparisons either confirm the
	 * key lies inside the predicted sixteenth — narrowing {@code [lo, hi]} to it and
	 * consuming the nibble — or push the interval to the side the comparisons indicate,
	 * which keeps the result exact for any (also non-uniform) key distribution. Once the
	 * interval gets small the search finishes with binary steps.
	 * <p>
	 * Compared to interpolation search the probe positions depend only on the key bits
	 * and {@code [lo, hi]} (not on the stored values), so the early, cache-critical probes
	 * land on a small fixed set of array positions that stay warm in the L3 cache, while
	 * the interval shrinks by ~16x per step instead of 2x — roughly halving the probe count.
	 *
	 * @return the index of {@code key} if found, otherwise {@code -(insertionPoint) - 1}
	 *         (same contract as {@link Arrays#binarySearch(long[], int, int, long)}).
	 */
	private int radixSearchSmall(final long key, final int to) {
		int lo = 0;
		int hi = to - 1;
		int shift = (k << 1) - 4; // top nibble of the 2*k significant bits
		while (lo <= hi) {
			final int span = hi - lo;
			if (shift < 0 || span < RADIX_MIN_SPAN) {
				final int mid = (lo + hi) >>> 1;
				final long midVal = kmers[mid];
				if (midVal < key) lo = mid + 1;
				else if (midVal > key) hi = mid - 1;
				else return mid;
				continue;
			}
			final int digit = (int) ((key >>> shift) & 0xF);
			final int lower = lo + (int) (((long) span * digit) >>> 4); // start of the d-th sixteenth
			final long lowerVal = kmers[lower];
			if (lowerVal == key) {
				return lower;
			}
			if (lowerVal > key) {
				hi = lower - 1; // key is before the predicted sixteenth
				continue;
			}
			final int upper = lo + (int) (((long) span * (digit + 1)) >>> 4); // start of the next sixteenth
			final long upperVal = kmers[upper];
			if (upperVal == key) {
				return upper;
			}
			if (upperVal < key) {
				lo = upper + 1; // key is beyond the predicted sixteenth
			} else {
				lo = lower + 1; // key is inside the predicted sixteenth: narrow and consume the nibble
				hi = upper - 1;
				shift -= 4;
			}
		}
		return -(lo + 1);
	}

	/**
	 * Radix-guided search over the large ({@code long[][]} BigArray) sorted k-mer store.
	 * The BigArray segment size is {@code 2^27}, so an index is split into segment
	 * ({@code pos >>> 27}) and offset ({@code pos & (2^27 - 1)}) via bit operations.
	 *
	 * @see #radixSearchSmall(long, int)
	 */
	private long radixSearchLarge(final long key, final long to) {
		long lo = 0;
		long hi = to - 1;
		int shift = (k << 1) - 4; // top nibble of the 2*k significant bits
		while (lo <= hi) {
			final long span = hi - lo;
			if (shift < 0 || span < RADIX_MIN_SPAN) {
				final long mid = (lo + hi) >>> 1;
				final long midVal = largeKmers[(int) (mid >>> 27)][(int) (mid & 134217727)];
				if (midVal < key) lo = mid + 1;
				else if (midVal > key) hi = mid - 1;
				else return mid;
				continue;
			}
			final int digit = (int) ((key >>> shift) & 0xF);
			final long lower = lo + ((span * digit) >>> 4); // start of the d-th sixteenth
			final long lowerVal = largeKmers[(int) (lower >>> 27)][(int) (lower & 134217727)];
			if (lowerVal == key) {
				return lower;
			}
			if (lowerVal > key) {
				hi = lower - 1; // key is before the predicted sixteenth
				continue;
			}
			final long upper = lo + ((span * (digit + 1)) >>> 4); // start of the next sixteenth
			final long upperVal = largeKmers[(int) (upper >>> 27)][(int) (upper & 134217727)];
			if (upperVal == key) {
				return upper;
			}
			if (upperVal < key) {
				lo = upper + 1; // key is beyond the predicted sixteenth
			} else {
				lo = lower + 1; // key is inside the predicted sixteenth: narrow and consume the nibble
				hi = upper - 1;
				shift -= 4;
			}
		}
		return -(lo + 1);
	}

	// Made final for potential (automated) inlining by JVM
	public final V getLong(final long kmer, final long[] posStore) {
		if (filter != null && useFilter && !filter.containsLong(kmer)) {
			return null;
		}
		short index;
		long pos;
		if (sorted) {
			if (largeKmers != null) {
				pos = radixSearch ? radixSearchLarge(kmer, entries)
						: LongBigArrays.binarySearch(largeKmers, 0, entries, kmer);
				if (pos < 0) {
					return null;
				}
				index = BigArrays.get(largeValueIndexes, pos);
			} else {
				pos = radixSearch ? radixSearchSmall(kmer, (int) entries)
						: Arrays.binarySearch(kmers, 0, (int) entries, kmer);
				if (pos < 0) {
					return null;
				}
				index = valueIndexes[(int) pos];
			}
		} else {
			if (largeKmers != null) {
				pos = -1;
				for (long i = 0; i < entries; i++) {
					if (BigArrays.get(largeKmers, i) == kmer) {
						pos = i;
						break;
					}
				}
				if (pos < 0) {
					return null;
				}
				index = BigArrays.get(largeValueIndexes, pos);
			} else {
				pos = -1;
				for (int i = 0; i < entries; i++) {
					if (kmers[i] == kmer) {
						pos = i;
						break;
					}
				}
				if (pos < 0) {
					return null;
				}
				index = valueIndexes[(int) pos];
			}
		}
		if (posStore != null) {
			posStore[0] = pos;
		}
		return indexMap[index - Short.MIN_VALUE];
	}

	// Made final for potential (automated) inlining by JVM
	public final V getLongInlined(final long kmer, final long[] posStore) {
		if (filter != null && useFilter && !filter.containsLong(kmer)) {
			return null;
		}
		short index;
		if (sorted) {
			if (largeKmers != null) {
				long pos;
				if (radixSearch) {
					pos = radixSearchLarge(kmer, entries);
				} else {
					pos = -1;
					boolean finished = false;
					long from = 0;
					long to = entries;
					long midVal;
					to--;
					while (from <= to) {
						final long mid = (from + to) >>> 1;
						midVal = largeKmers[(int) (mid >>> 27)][(int) (mid & 134217727)];
						if (midVal < kmer) from = mid + 1;
						else if (midVal > kmer) to = mid - 1;
						else {
							pos = mid;
							finished = true;
							break;
						}
					}
					if (!finished) {
						pos = -(from + 1);
					}
				}
				if (posStore != null) {
					posStore[0] = pos;
				}
				if (pos < 0) {
					return null;
				}
				index = largeValueIndexes[(int) (pos >>> 27)][(int) (pos & 134217727)];
			} else {
				int pos;
				if (radixSearch) {
					pos = radixSearchSmall(kmer, (int) entries);
				} else {
					pos = -1;
					int low = 0;
					int high = (int) entries - 1;
					boolean finished = false;

					while (low <= high) {
						int mid = (low + high) >>> 1;
						long midVal = kmers[mid];

						if (midVal < kmer)
							low = mid + 1;
						else if (midVal > kmer)
							high = mid - 1;
						else {
							pos = mid;
							finished = true;
							break;// key found
						}
					}
					if (!finished) {
						pos = -(low + 1);// key not found.
					}
				}
				if (posStore != null) {
					posStore[0] = pos;
				}
				if (pos < 0) {
					return null;
				}
				index = valueIndexes[pos];
			}
		} else {
			if (largeKmers != null) {
				long pos = -1;
				for (long i = 0; i < entries; i++) {
					if (BigArrays.get(largeKmers, i) == kmer) {
						pos = i;
						break;
					}
				}
				if (posStore != null) {
					posStore[0] = pos;
				}
				if (pos < 0) {
					return null;
				}
				index = BigArrays.get(largeValueIndexes, pos);
			} else {
				int pos = -1;
				for (int i = 0; i < entries; i++) {
					if (kmers[i] == kmer) {
						pos = i;
						break;
					}
				}
				if (posStore != null) {
					posStore[0] = pos;
				}
				if (pos < 0) {
					return null;
				}
				index = valueIndexes[pos];
			}
		}
		return indexMap[index - Short.MIN_VALUE];
	}

	public int indexAtPosition(long pos) {
		return (largeKmers != null ? BigArrays.get(largeValueIndexes, pos) : valueIndexes[(int) pos]) - Short.MIN_VALUE;
	}

	@Override
	public void optimize() {
		if (sorted) {
			return;
		}
		if (largeKmers != null) {
			BigArrays.quickSort(0, entries, new LongComparator() {
				@Override
				public int compare(long k1, long k2) {
					long kmer1 = BigArrays.get(largeKmers, k1);
					long kmer2 = BigArrays.get(largeKmers, k2);
					return Long.compare(kmer1, kmer2);
				}
			}, new BigSwapper() {
				@Override
				public void swap(long a, long b) {
					long kmerA = BigArrays.get(largeKmers, a);
					long kmerB = BigArrays.get(largeKmers, b);
					BigArrays.set(largeKmers, b, kmerA);
					BigArrays.set(largeKmers, a, kmerB);
					short indexA = BigArrays.get(largeValueIndexes, a);
					short indexB = BigArrays.get(largeValueIndexes, b);
					BigArrays.set(largeValueIndexes, b, indexA);
					BigArrays.set(largeValueIndexes, a, indexB);
				}
			});
		} else {
			BigArrays.quickSort(0, entries, new LongComparator() {
				@Override
				public int compare(long k1, long k2) {
					long kmer1 = kmers[(int) k1];
					long kmer2 = kmers[(int) k2];
					return Long.compare(kmer1, kmer2);
				}
			}, new BigSwapper() {
				@Override
				public void swap(long a, long b) {
					long kmerA = kmers[(int) a];
					long kmerB = kmers[(int) b];
					kmers[(int) b] = kmerA;
					kmers[(int) a] = kmerB;
					short indexA = valueIndexes[(int) a];
					short indexB = valueIndexes[(int) b];
					valueIndexes[(int) b] = indexA;
					valueIndexes[(int) a] = indexB;
				}
			});
		}
		sorted = true;
		// Rework the bloom filter...
		boolean xor = filter instanceof XORKMerBloomFilter;
		filter = null; // Set to null already for garbage collector.
		if (optimizedFpp < 1) {
			if (optimizedFpp == BlockedKMerBloomFilter.DEFAULT_FPP) {
				filter = new BlockedKMerBloomFilter();
			}
			else {
				filter = xor ? new XORKMerBloomFilter(optimizedFpp) : new MurmurKMerBloomFilter(optimizedFpp);
			}
			// Now make a filter with High fpp.
			filter.ensureExpectedSize(entries, false);
			long kmer;
			for (long i = 0; i < entries; i++) {
				if (largeKmers != null) {
					kmer = BigArrays.get(largeKmers, i);
				} else {
					kmer = kmers[(int) i];
				}
				filter.putLong(kmer);
			}
		}
	}

	@Override
	public boolean isOptimized() {
		return sorted;
	}

	@Override
	public void visit(KMerStoreVisitor<V> visitor) {
		V value;
		long kmer;
		short sindex;
		for (long i = 0; i < entries; i++) {
			if (largeKmers != null) {
				kmer = BigArrays.get(largeKmers, i);
				sindex = BigArrays.get(largeValueIndexes, i);
			} else {
				kmer = kmers[(int) i];
				sindex = valueIndexes[(int) i];
			}
			value = indexMap[sindex - Short.MIN_VALUE];
			visitor.nextValue(this, kmer, value);
		}
	}

	public void visit(KMerSortedArrayVisitor<V> visitor) {
		long kmer;
		int sindex;
		for (long i = 0; i < entries; i++) {
			if (largeKmers != null) {
				kmer = BigArrays.get(largeKmers, i);
				sindex = BigArrays.get(largeValueIndexes, i);
			} else {
				kmer = kmers[(int) i];
				sindex = valueIndexes[(int) i];
			}
			visitor.nextValue(this, kmer, sindex  - Short.MIN_VALUE, i);
		}
	}

	public interface KMerSortedArrayVisitor<V extends Serializable> {
		public void nextValue(KMerSortedArray<V> trie, long kmer, int index, long i);
	}

	public interface UpdateValueProvider<V extends Serializable> {
		public V getUpdateValue(V oldValue);
	}

	public interface ValueConverter<V, W extends Serializable> {
		public W convertValue(V value);
	}
}
