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
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;

import org.metagene.genestrip.bloom.MurmurCGATBloomFilter;
import org.metagene.genestrip.util.CGAT;
import org.metagene.genestrip.util.CGATRingBuffer;

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

	public static byte EXCLUDED_KMER_VIA_COUNT = Byte.MIN_VALUE;

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

	protected long size;

	private long entries;
	private boolean sorted;
	private boolean initSize;
	private short nextValueIndex;
	private transient MurmurCGATBloomFilter filter;
	private boolean useFilter;
	private Object2LongMap<V> kmerPersTaxid;

	private transient long kmersMoved;

	public KMerSortedArray(int k, double fpp, List<V> initialValues, boolean enforceLarge) {
		this(k, fpp, initialValues, enforceLarge, new MurmurCGATBloomFilter(k, fpp));
	}

	@SuppressWarnings("unchecked")
	protected KMerSortedArray(int k, double fpp, List<V> initialValues, boolean enforceLarge,
			MurmurCGATBloomFilter filter) {
		this.k = k;
		int s = initialValues == null ? 0 : initialValues.size();
		indexMap = (V[]) new Serializable[MAX_VALUES];
		valueMap = new Object2ShortOpenHashMap<V>(s);
		// It is VERY important to start with one here because of a bug in the fastutil
		// library.
		// It seems to have to do with storing 0 as a key and deserializing
		// Short2ObjectHashMap in this case.
		// I did not want to dive into this - so this is a simple workaround:
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
		valueMap = new Object2ShortOpenHashMap<V>(org.valueMap.size());

		if (org.kmerPersTaxid != null) {
			kmerPersTaxid = new Object2LongOpenHashMap<V>();
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

	public MurmurCGATBloomFilter getFilter() {
		return filter;
	}

	public void setFilter(MurmurCGATBloomFilter filter) {
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

		Object2LongMap<V> map = new Object2LongOpenHashMap<V>();
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
	public boolean put(CGATRingBuffer buffer, V value, boolean reverse) {
		throw new UnsupportedOperationException("Deprecated and almost removed.");
	}

	@Override
	public boolean put(byte[] nseq, int start, V value, boolean reverse) {
		long kmer = reverse ? CGAT.kMerToLongReverse(nseq, start, k, null)
				: CGAT.kMerToLongStraight(nseq, start, k, null);
		return putLong(kmer, value);
	}

	public boolean isFull() {
		return entries == size;
	}

	/**
	 * @return true if put succeeded, false if (probably) some value is already
	 *         stored under that kmer.
	 */
	public boolean putLong(long kmer, V value) {
		if (value == null) {
			throw new NullPointerException("null is not allowed as a value.");
		}
		if (filter.containsLong(kmer)) {
			// Fail fast - we could check if the kmer is indeed stored, but it's way too
			// slow because
			// of linear search in the kmer array...
			return false;
		}
		sorted = false;
		if (entries == size) {
			throw new IllegalStateException("Capacity exceeded.");
		}
		short sindex = getAddValueIndex(value);
		if (largeKmers != null) {
			BigArrays.set(largeKmers, entries, kmer);
			BigArrays.set(largeValueIndexes, entries, sindex);
		} else {
			kmers[(int) entries] = kmer;
			valueIndexes[(int) entries] = sindex;
		}
		entries++;
		filter.putLong(kmer);
		return true;
	}

	protected short getAddValueIndex(V value) {
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
	public V get(CGATRingBuffer buffer, boolean reverse) {
		throw new UnsupportedOperationException("Deprecated and almost removed.");
	}

	@Override
	public V get(byte[] nseq, int start, boolean reverse) {
		long kmer = reverse ? CGAT.kMerToLongReverse(nseq, start, k, null)
				: CGAT.kMerToLongStraight(nseq, start, k, null);
		return getLong(kmer, null);
	}

	public V get(byte[] nseq, int start, boolean reverse, long[] indexStore) {
		long kmer = reverse ? CGAT.kMerToLongReverse(nseq, start, k, null)
				: CGAT.kMerToLongStraight(nseq, start, k, null);
		return getLong(kmer, indexStore);
	}

	public boolean update(long kmer, UpdateValueProvider<V> provider) {
		if (!sorted) {
			throw new IllegalStateException("Updated only works when optimized.");
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
		synchronized (provider.getSynchronizationObject(pos)) {
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
				if (largeKmers != null) {
					BigArrays.set(largeValueIndexes, pos, index);
				} else {
					valueIndexes[(int) pos] = index;
				}
				return true;
			}
			return false;
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
	// Also inlined the sorted, largeKmers != null part for speed.
	public final V getLongInlined(final long kmer, final long[] posStore) {
		if (filter != null && useFilter && !filter.containsLong(kmer)) {
			return null;
		}
		short index;
		long pos = 0;
		if (sorted) {
			if (largeKmers != null) {
				boolean finished = false;
				long from = 0;
				long to = entries;
				long midVal;
				to--;
				while (from <= to) {
					final long mid = (from + to) >>> 1;
					midVal = largeKmers[(int) (mid >>> BigArrays.SEGMENT_SHIFT)][(int) (mid & BigArrays.SEGMENT_MASK)];
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
				if (pos < 0) {
					return null;
				}
				index = largeValueIndexes[(int) (pos >>> BigArrays.SEGMENT_SHIFT)][(int) (pos & BigArrays.SEGMENT_MASK)];
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

		public Object getSynchronizationObject(long position);
	}

	public interface ValueConverter<V, W extends Serializable> {
		public W convertValue(V value);
	}
}
