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

import org.metagene.genestrip.bloom.AbstractKMerBloomFilter;
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
	public static final int MAX_VALUES = Short.MAX_VALUE;

	private static final long serialVersionUID = 1L;

	public static long MAX_SMALL_CAPACITY = Integer.MAX_VALUE - 8;

	private long[] kmers;
	private short[] valueIndexes;

	private long[][] largeKmers;
	private short[][] largeValueIndexes;

	private final int k;
	private final Object2ShortMap<V> valueMap;
	// For efficiency and synchronization problems this was turned into an array of
	// lenth Short.MAX_VALUE. Only 250kb in size...
	private final V[] indexMap;
	private final boolean enforceLarge;

	// Just for optimizing synchronization during updates
	private transient Object[] syncs;

	protected long size;

	private long entries;
	private boolean sorted;
	private boolean initSize;
	private short nextValueIndex;
	private transient AbstractKMerBloomFilter filter;
	private boolean useFilter;
	private Object2LongMap<V> kmerPersTaxid;

	private transient long kmersMoved;

	public KMerSortedArray(int k, double fpp, List<V> initialValues, boolean enforceLarge, boolean xor) {
		this(k, initialValues, enforceLarge, xor ? new XORKMerBloomFilter(k, fpp) : new MurmurKMerBloomFilter(k, fpp));
	}

	@SuppressWarnings("unchecked")
	protected KMerSortedArray(int k, List<V> initialValues, boolean enforceLarge,
			AbstractKMerBloomFilter filter) {
		this.k = k;
		int s = initialValues == null ? 0 : initialValues.size();
		indexMap = (V[]) new Serializable[MAX_VALUES];
		valueMap = new Object2ShortOpenHashMap<>(s);
		initSyncs();
		nextValueIndex = 0;
		this.enforceLarge = enforceLarge;
		if (initialValues != null) {
			for (V v : initialValues) {
				getAddValueIndex(v);
			}
		}
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

		size = org.size;

		entries = org.entries;
		sorted = org.sorted;
		initSize = org.initSize;
		nextValueIndex = org.nextValueIndex;
		filter = org.filter;
		useFilter = org.useFilter;

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
			private int i = 0;

			@Override
			public boolean hasNext() {
				return i < nextValueIndex;
			}

			@Override
			public V next() {
				return indexMap[i++];
			}
		};
	}

	public AbstractKMerBloomFilter getFilter() {
		return filter;
	}

	public void setFilter(AbstractKMerBloomFilter filter) {
		this.filter = filter;
	}

	public void setUseFilter(boolean useFilter) {
		this.useFilter = useFilter;
	}

	public boolean isUseFilter() {
		return useFilter;
	}

	public Object2LongMap<V> getNKmersPerTaxid() {
		long[] countArray = new long[getNValues()];
		if (largeKmers != null) {
			for (long i = 0; i < entries; i++) {
				countArray[BigArrays.get(largeValueIndexes, i)]++;
			}
		} else {
			for (int i = 0; i < entries; i++) {
				countArray[valueIndexes[i]]++;
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
		return nextValueIndex;
	}

	public V getValueForIndex(short index) {
		return indexMap[index];
	}

	public short getIndexForValue(V value) {
		return valueMap.getOrDefault(value, (short) -1);
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
		short sindex;
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

	public short getAddValueIndex(V value) {
		short sindex = valueMap.getOrDefault(value, (short) -1);
		if (sindex == -1) {
			if (nextValueIndex == MAX_VALUES) {
				throw new IllegalStateException("Too many different values - only " + MAX_VALUES + " are possible.");
			}
			valueMap.put(value, nextValueIndex);
			indexMap[nextValueIndex] = value;
			sindex = nextValueIndex;
			nextValueIndex++;
		}
		return sindex;
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
			pos = LongBigArrays.binarySearch(largeKmers, 0, entries, kmer);
			if (pos < 0) {
				return false;
			}
		} else {
			pos = Arrays.binarySearch(kmers, 0, (int) entries, kmer);
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
			short index;
			if (largeKmers != null) {
				index = BigArrays.get(largeValueIndexes, pos);
			} else {
				index = valueIndexes[(int) pos];
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
				setIndexAtPosition(kmer, index);
				return true;
			}
			return false;
		}
	}

	public void setIndexAtPosition(long pos, short index) {
		if (largeKmers != null) {
			BigArrays.set(largeValueIndexes, pos, index);
		} else {
			valueIndexes[(int) pos] = index;
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

	// Made final for potential (automated) inlining by JVM
	public final V getLong(final long kmer, final long[] posStore) {
		if (filter != null && useFilter && !filter.containsLong(kmer)) {
			return null;
		}
		short index;
		long pos;
		if (sorted) {
			if (largeKmers != null) {
				pos = LongBigArrays.binarySearch(largeKmers, 0, entries, kmer);
				if (pos < 0) {
					return null;
				}
				index = BigArrays.get(largeValueIndexes, pos);
			} else {
				pos = Arrays.binarySearch(kmers, 0, (int) entries, kmer);
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
		return indexMap[index];
	}

	// Made final for potential (automated) inlining by JVM
	public final V getLongInlined(final long kmer, final long[] posStore) {
		if (filter != null && useFilter && !filter.containsLong(kmer)) {
			return null;
		}
		short index;
		if (sorted) {
			if (largeKmers != null) {
				long pos = -1;
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
				if (posStore != null) {
					posStore[0] = pos;
				}
				if (pos < 0) {
					return null;
				}
				index = largeValueIndexes[(int) (pos >>> 27)][(int) (pos & 134217727)];
			} else {
				int pos = -1;
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
		return indexMap[index];
	}

	public short indexAtPosition(long pos) {
		return largeKmers != null ? BigArrays.get(largeValueIndexes, pos) : valueIndexes[(int) pos];
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
			value = indexMap[sindex];
			visitor.nextValue(this, kmer, value);
		}
	}

	public void visit(KMerSortedArrayVisitor<V> visitor) {
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
			visitor.nextValue(this, kmer, sindex, i);
		}
	}

	public interface KMerSortedArrayVisitor<V extends Serializable> {
		public void nextValue(KMerSortedArray<V> trie, long kmer, short index, long i);
	}

	public interface UpdateValueProvider<V extends Serializable> {
		public V getUpdateValue(V oldValue);
	}

	public interface ValueConverter<V, W extends Serializable> {
		public W convertValue(V value);
	}
}
