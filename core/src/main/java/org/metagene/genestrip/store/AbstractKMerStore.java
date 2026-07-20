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
import java.io.ObjectOutputStream;
import java.io.Serializable;
import java.util.Iterator;
import java.util.List;

import org.metagene.genestrip.bloom.BlockedKMerBloomFilter;
import org.metagene.genestrip.bloom.KMerProbFilter;
import org.metagene.genestrip.bloom.MurmurKMerBloomFilter;
import org.metagene.genestrip.bloom.XORKMerBloomFilter;

import it.unimi.dsi.fastutil.objects.Object2LongMap;
import it.unimi.dsi.fastutil.objects.Object2LongOpenHashMap;
import it.unimi.dsi.fastutil.objects.Object2IntMap;
import it.unimi.dsi.fastutil.objects.Object2IntOpenHashMap;

/**
 * Common base for the {@link KMerStore} implementations {@link KMerSortedArray} and
 * {@link RadixKMerStore}. It owns everything that does not depend on the concrete k-mer
 * storage layout:
 * <ul>
 * <li>the value &lt;-&gt; store-index mapping ({@link #valueMap} / {@link #indexMap} and the
 *     {@code getAddValueIndex} / {@code getValueForIndex} / {@code getIndexForValue} /
 *     {@code getNValues} / {@code getValues} accessors),</li>
 * <li>the basic counters ({@code k}, {@code entries}, {@code size}, {@code sorted}) and their
 *     getters,</li>
 * <li>the probabilistic pre-filter ({@code getFilter} / {@code setFilter} / {@code setUseFilter} /
 *     {@code isUseFilter}) plus the {@link #createOptimizedFilter()} helper used when sorting,</li>
 * <li>the per-value k-mer count statistics ({@code getNKmersPerTaxid} /
 *     {@code invalidateNKmersPerTaxid}, cached and baked into the serialized form on
 *     {@code writeObject}).</li>
 * </ul>
 * The storage-specific operations ({@code putLong}, {@code getLong}, {@code update},
 * {@code optimize}, {@code visit}, {@code convertValues}) remain abstract.
 *
 * @param <V> the value type mapped to each k-mer
 */
public abstract class AbstractKMerStore<V extends Serializable> implements KMerStore<V> {
	private static final long serialVersionUID = 1L;

	/** The k-mer length; always in {@code [1, 31]}. */
	protected final int k;

	// Maximum number of distinct values this store supports. Defined per subclass (see their
	// MAX_VALUES constants) since it depends on how the value index is stored (short array vs.
	// packed bits) and passed to this base constructor.
	/** The value-index capacity (maximum number of distinct values). */
	protected final int maxValues;
	// Maps a value to its 0-based store index.
	/** Maps each value to its 0-based store index. */
	protected final Object2IntMap<V> valueMap;
	// Indexed 0-based by value index; length maxValues.
	/** Maps each store index back to its value; length {@code maxValues}. */
	protected final V[] indexMap;
	/** The next store index to assign to a newly registered value. */
	protected int nextValueIndex;

	/** The total capacity of the store. */
	protected long size;
	/** The number of entries currently stored. */
	protected long entries;
	/** Whether the store has been sorted/optimized. */
	protected boolean sorted;

	/** The target false-positive probability of the post-optimize filter. */
	protected final double optimizedFpp;
	// Transient: the filter is (re)built on optimize() and stored separately by Database.
	/** The optional probabilistic pre-filter, or {@code null} if none. */
	protected transient KMerProbFilter filter;
	/** Whether the pre-filter is currently in use. */
	protected boolean useFilter;

	/** Per-value k-mer count statistics, or {@code null} if not tracked. */
	protected Object2LongMap<V> kmerPersTaxid;

	/**
	 * Base constructor: sets {@code k} (which must be in {@code [1, 31]}), the value-index capacity,
	 * the optional fill-time pre-filter and the post-optimize filter fpp, and registers any initial
	 * values so they receive stable store indices.
	 *
	 * @param k             the k-mer length; must be in {@code [1, 31]}
	 * @param maxValues     the value-index capacity (maximum number of distinct values)
	 * @param initialValues values to register up front for stable store indices, or {@code null}
	 * @param filter        the optional fill-time pre-filter, or {@code null} for none
	 * @param optimizedFpp  the target false-positive probability of the post-optimize filter
	 */
	@SuppressWarnings("unchecked")
	protected AbstractKMerStore(int k, int maxValues, List<V> initialValues, KMerProbFilter filter, double optimizedFpp) {
		// k must be <= 31: at k=32 a k-mer fills all 64 bits, so the all-T k-mer collides with the
		// -1L "invalid k-mer" sentinel used throughout CGAT/the matcher (see GSConfigKey.KMER_SIZE).
		if (k < 1 || k > 31) {
			throw new IllegalArgumentException("k must be in [1, 31], got " + k);
		}
		this.k = k;
		this.maxValues = maxValues;
		int s = initialValues == null ? 0 : initialValues.size();
		indexMap = (V[]) new Serializable[maxValues];
		valueMap = new Object2IntOpenHashMap<>(s);
		nextValueIndex = 0;
		if (initialValues != null) {
			for (V v : initialValues) {
				getAddValueIndex(v);
			}
		}
		this.optimizedFpp = optimizedFpp;
		this.filter = filter;
		this.useFilter = filter != null;
		this.kmerPersTaxid = null;
	}

	/**
	 * Copy constructor shared by the subclasses' {@code convertValues} support: copies the common
	 * counters and rebuilds the value maps with converted values. The concrete k-mer storage is
	 * copied (typically shared) by the subclass constructor.
	 *
	 * @param <W>       the value type of the store being converted from
	 * @param org       the store whose common state is copied
	 * @param converter converts each value of {@code org} to a value of this store
	 */
	@SuppressWarnings("unchecked")
	protected <W extends Serializable> AbstractKMerStore(AbstractKMerStore<W> org,
														 KMerStore.ValueConverter<W, V> converter) {
		k = org.k;
		maxValues = org.maxValues;
		optimizedFpp = org.optimizedFpp;
		size = org.size;
		entries = org.entries;
		sorted = org.sorted;
		nextValueIndex = org.nextValueIndex;
		filter = org.filter;
		useFilter = org.useFilter;

		indexMap = (V[]) new Serializable[maxValues];
		valueMap = new Object2IntOpenHashMap<>(org.valueMap.size());
		if (org.kmerPersTaxid != null) {
			kmerPersTaxid = new Object2LongOpenHashMap<>();
		}
		int nValues = org.getNValues();
		for (int i = 0; i < nValues; i++) {
			W orgValue = org.indexMap[i];
			V value = converter.convertValue(orgValue);
			indexMap[i] = value;
			valueMap.put(value, i);
			if (kmerPersTaxid != null) {
				kmerPersTaxid.put(value, org.kmerPersTaxid.getLong(orgValue));
			}
		}
	}

	@Override
	public int getK() {
		return k;
	}

	@Override
	public int getMaxValues() {
		return maxValues;
	}

	@Override
	public long getEntries() {
		return entries;
	}

	@Override
	public long getSize() {
		return size;
	}

	@Override
	public boolean isFull() {
		return entries == size;
	}

	@Override
	public boolean isOptimized() {
		return sorted;
	}

	// --- Probabilistic pre-filter ---------------------------------------------

	@Override
	public KMerProbFilter getFilter() {
		return filter;
	}

	@Override
	public void setFilter(KMerProbFilter filter) {
		this.filter = filter;
	}

	@Override
	public void setUseFilter(boolean useFilter) {
		this.useFilter = useFilter;
	}

	@Override
	public boolean isUseFilter() {
		return useFilter;
	}

	/**
	 * Builds the post-{@code optimize()} pre-filter sized for {@link #entries}, choosing the same
	 * kind (blocked / XOR / Murmur) as the current fill-time filter. Returns {@code null} when
	 * {@code optimizedFpp >= 1} (no filter). The caller is expected to assign the result to
	 * {@link #filter} and then populate it with the stored k-mers.
	 *
	 * @return the newly created pre-filter, or {@code null} if {@code optimizedFpp >= 1}
	 */
	protected KMerProbFilter createOptimizedFilter() {
		if (optimizedFpp >= 1) {
			return null;
		}
		KMerProbFilter f;
		if (optimizedFpp == BlockedKMerBloomFilter.DEFAULT_FPP) {
			f = new BlockedKMerBloomFilter(entries);
		} else {
			boolean xor = filter instanceof XORKMerBloomFilter;
			f = xor ? new XORKMerBloomFilter(optimizedFpp, entries) : new MurmurKMerBloomFilter(optimizedFpp, entries);
		}
		return f;
	}

	// --- Value <-> store-index mapping ----------------------------------------

	@Override
	public int getNValues() {
		return nextValueIndex;
	}

	@Override
	public V getValueForIndex(int index) {
		return indexMap[index];
	}

	@Override
	public int getIndexForValue(V value) {
		return valueMap.getOrDefault(value, -1);
	}

	@Override
	public int getAddValueIndex(V value) {
		int index = valueMap.getOrDefault(value, -1);
		if (index < 0) {
			if (nextValueIndex >= maxValues) {
				throw new IllegalStateException("Too many different values - only " + maxValues + " are possible.");
			}
			index = nextValueIndex++;
			valueMap.put(value, index);
			indexMap[index] = value;
		}
		return index;
	}

	/**
	 * Resolves an <em>already registered</em> value to its store index with a lock-free read, for the
	 * concurrent fill and update paths. It never mutates {@code valueMap} (a plain, non-thread-safe hash
	 * map), so it needs no synchronization: both phases run only after every value has been registered
	 * up front (single-threaded — the fill via the store's {@code initialValues}, the update via {@link
	 * Database#initStoreIndices()}), so no concurrent writer can race the reads. Encountering an
	 * unregistered value would mean that invariant was violated; rather than silently mutate the map
	 * from many threads (which is unsafe) it fails fast.
	 *
	 * @param value the value to resolve; must already be registered
	 * @return the store index of the value
	 * @throws IllegalStateException if the value was not registered up front
	 */
	protected final int getRegisteredValueIndex(V value) {
		int index = getIndexForValue(value);
		if (index < 0) {
			throw new IllegalStateException("Encountered a value not registered before the concurrent phase: "
					+ value + ". All values must be registered up front (fill: the store's initialValues;"
					+ " update: Database.initStoreIndices()).");
		}
		return index;
	}

	@Override
	public Iterator<V> getValues() {
		return new Iterator<V>() {
			private int i = 0;
			private final int n = getNValues();

			@Override
			public boolean hasNext() {
				return i < n;
			}

			@Override
			public V next() {
				return indexMap[i++];
			}
		};
	}

	// --- Statistics -----------------------------------------------------------

	@Override
	public Object2LongMap<V> getNKmersPerTaxid() {
		// Lazily cached; the value-mutating operations (optimize() and setIndexAtPosition()) invalidate
		// the cache, so a repeated read is free and a stale count is never returned.
		if (this.kmerPersTaxid == null) {
			this.kmerPersTaxid = computeNKmersPerTaxid();
		}
		return this.kmerPersTaxid;
	}

	/**
	 * Counts, by visiting every entry, how many stored k-mers map to each value (the {@code null} key
	 * maps to the total number of entries). Used to fill the {@link #kmerPersTaxid} cache.
	 *
	 * @return a freshly computed map from each value to its number of stored k-mers.
	 */
	private Object2LongMap<V> computeNKmersPerTaxid() {
		final long[] countArray = new long[getNValues()];
		visit(new KMerStore.IndexedKMerStoreVisitor<V>() {
			@Override
			public void nextValue(KMerStore<V> store, long kmer, int index, long pos) {
				countArray[index]++;
			}
		});

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

	@Override
	public final void invalidateNKmersPerTaxid() {
		// Guarded so that the hot, concurrent setIndexAtPosition() path only writes the field once per
		// invalidation rather than on every call, keeping the shared field out of cache-line ping-pong.
		if (this.kmerPersTaxid != null) {
			this.kmerPersTaxid = null;
		}
	}

	/**
	 * Ensures the per-value k-mer counts are baked in before the store is serialized, so a loaded store
	 * returns them without revisiting every entry. Filling the cache if absent suffices: value
	 * reassignments invalidate it (automatically via {@link #setIndexAtPosition(long, int)}, or
	 * explicitly via {@link #invalidateNKmersPerTaxid()} after a bulk update), so it is never stale here.
	 *
	 * @param out the stream to write to
	 * @throws IOException if writing fails
	 */
	private void writeObject(ObjectOutputStream out) throws IOException {
		getNKmersPerTaxid();
		out.defaultWriteObject();
	}
}
