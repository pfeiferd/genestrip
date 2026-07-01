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
 * <li>the optional probabilistic pre-filter ({@link TunableKMerStore}: {@code getFilter} /
 *     {@code setFilter} / {@code setUseFilter} / {@code isUseFilter}) plus the
 *     {@link #createOptimizedFilter()} helper used when sorting,</li>
 * <li>the per-value k-mer count statistics ({@code fix} / {@code getFixedNKmersPerTaxid} /
 *     {@code getNKmersPerTaxid}), and</li>
 * <li>the lock array used to serialize concurrent updates of the same entry.</li>
 * </ul>
 * The storage-specific operations ({@code initSize}, {@code putLong}, {@code getLong},
 * {@code update}, {@code optimize}, {@code visit}, {@code convertValues}) remain abstract.
 */
public abstract class AbstractKMerStore<V extends Serializable> implements TunableKMerStore<V> {
	private static final long serialVersionUID = 1L;

	protected final int k;

	// Maximum number of distinct values this store supports. Defined per subclass (see their
	// MAX_VALUES constants) since it depends on how the value index is stored (short array vs.
	// packed bits) and passed to this base constructor.
	protected final int maxValues;
	// Maps a value to its 0-based store index.
	protected final Object2IntMap<V> valueMap;
	// Indexed 0-based by value index; length maxValues.
	protected final V[] indexMap;
	protected int nextValueIndex;

	protected long size;
	protected long entries;
	protected boolean sorted;

	protected final double optimizedFpp;
	// Transient: the filter is (re)built on optimize() and stored separately by Database.
	protected transient KMerProbFilter filter;
	protected boolean useFilter;

	protected Object2LongMap<V> kmerPersTaxid;

	protected transient long kmersMoved;
	// Just for optimizing synchronization during updates.
	protected transient Object[] syncs;

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
		this.kmersMoved = 0;
		initSyncs();
	}

	/**
	 * Copy constructor shared by the subclasses' {@code convertValues} support: copies the common
	 * counters and rebuilds the value maps with converted values. The concrete k-mer storage is
	 * copied (typically shared) by the subclass constructor.
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
		initSyncs();
	}

	private void readObject(ObjectInputStream ois) throws ClassNotFoundException, IOException {
		ois.defaultReadObject();
		initSyncs();
		kmersMoved = 0;
	}

	private void initSyncs() {
		syncs = new Object[512];
		for (int i = 0; i < syncs.length; i++) {
			syncs[i] = new Object();
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

	@Override
	public long getKMersMoved() {
		return kmersMoved;
	}

	// --- Optional pre-filter (TunableKMerStore) -------------------------------

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
	 */
	protected KMerProbFilter createOptimizedFilter() {
		if (optimizedFpp >= 1) {
			return null;
		}
		KMerProbFilter f;
		if (optimizedFpp == BlockedKMerBloomFilter.DEFAULT_FPP) {
			f = new BlockedKMerBloomFilter();
		} else {
			boolean xor = filter instanceof XORKMerBloomFilter;
			f = xor ? new XORKMerBloomFilter(optimizedFpp) : new MurmurKMerBloomFilter(optimizedFpp);
		}
		f.ensureExpectedSize(entries, false);
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

	public Object2LongMap<V> getNKmersPerTaxid() {
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
	public void fix() {
		this.kmerPersTaxid = getNKmersPerTaxid();
	}

	@Override
	public Object2LongMap<V> getFixedNKmersPerTaxid() {
		if (this.kmerPersTaxid == null) {
			fix();
		}
		return this.kmerPersTaxid;
	}
}
