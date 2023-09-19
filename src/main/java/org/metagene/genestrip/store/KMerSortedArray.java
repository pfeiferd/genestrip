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
import it.unimi.dsi.fastutil.shorts.Short2ObjectMap;
import it.unimi.dsi.fastutil.shorts.Short2ObjectOpenHashMap;

public class KMerSortedArray<V extends Serializable> implements KMerStore<V> {
	public static final int MAX_VALUES = Short.MAX_VALUE;

	private static final long serialVersionUID = 1L;

	public static long MAX_SMALL_CAPACITY = Integer.MAX_VALUE - 8;

	public static byte EXCLUDED_KMER_VIA_COUNT = Byte.MIN_VALUE;

	private long[] kmers;
	private byte[] counts;
	private short[] valueIndexes;

	private long[][] largeKmers;
	private byte[][] largeCounts;
	private short[][] largeValueIndexes;

	private final int k;
	private final Object2ShortMap<V> valueMap;
	private final Short2ObjectMap<V> indexMap;
	private final boolean enforceLarge;
	private final boolean withCounts;

	protected long size;

	private long entries;
	private boolean sorted;
	private boolean countsChanged;
	private boolean initSize;
	private short nextValueIndex;
	private MurmurCGATBloomFilter filter;

	public KMerSortedArray(int k, double fpp, List<V> initialValues, boolean enforceLarge, boolean withCounts,
			MurmurCGATBloomFilter filter) {
		this.k = k;
		this.withCounts = withCounts;
		int s = initialValues == null ? 0 : initialValues.size();
		indexMap = new Short2ObjectOpenHashMap<V>(s);
		valueMap = new Object2ShortOpenHashMap<V>(s);
		nextValueIndex = 0;
		this.enforceLarge = enforceLarge;
		if (initialValues != null) {
			for (V v : initialValues) {
				getAddValueIndex(v);
			}
		}
		this.filter = filter;
	}
	
	public long[] getStatsAsIndexArray() {		
		long[] counts = new long[nextValueIndex];
		visit(new KMerSortedArrayVisitor<V>() {
			@Override
			public void nextValue(KMerSortedArray<V> trie, long kmer, short index, long i) {
				counts[index]++;
			}
		});
		return counts;
	}
	
	public Object2LongMap<V> getStats() {		
		long[] counts = getStatsAsIndexArray();
		
		Object2LongMap<V> res = new Object2LongOpenHashMap<V>(counts.length);
		for (short i = 0; i < counts.length; i++) {
			res.put(indexMap.get(i), counts[i]);
		}
		res.put(null, entries);
		
		return res;
	}

	public KMerSortedArray(int k, double fpp, List<V> initialValues, boolean enforceLarge, boolean withCounts) {
		this(k, fpp, initialValues, enforceLarge, withCounts, new MurmurCGATBloomFilter(k, fpp));
	}

	public int getNValues() {
		return nextValueIndex;
	}

	public V getValueForIndex(short index) {
		return indexMap.get(index);
	}

	public short getIndexForValue(V value) {
		return valueMap.getOrDefault(value, (short) -1);
	}

	@Override
	public void initSize(long size) {
		if (initSize) {
			throw new IllegalStateException("Cant initlialize size twice.");
		}
		if (size <= 0) {
			throw new IllegalArgumentException("Expected insertions must be > 0.");
		}
		initSize = true;
		filter.clear();
		filter.ensureExpectedSize(size, false);
		this.size = size;

		if (size > MAX_SMALL_CAPACITY || enforceLarge) {
			largeKmers = BigArrays.ensureCapacity(BigArrays.wrap(new long[0]), size);
			largeValueIndexes = BigArrays.ensureCapacity(BigArrays.wrap(new short[0]), size);
			if (withCounts) {
				largeCounts = BigArrays.ensureCapacity(BigArrays.wrap(new byte[0]), size);
			}
		} else {
			kmers = new long[(int) size];
			valueIndexes = new short[(int) size];
			if (withCounts) {
				counts = new byte[(int) size];
			}
		}
	}

	protected void clearArray() {
		if (largeKmers != null) {
			BigArrays.fill(largeKmers, 0);
			BigArrays.fill(largeValueIndexes, (short) 0);
			if (withCounts) {
				BigArrays.fill(largeCounts, (byte) 0);
				countsChanged = false;
			}
		} else {
			Arrays.fill(kmers, 0);
			Arrays.fill(valueIndexes, (short) 0);
			if (withCounts) {
				Arrays.fill(counts, (byte) 0);
				countsChanged = false;
			}
		}
	}

	public Object2LongMap<V> getNKmersPerTaxid() {
		long[] countArray = new long[nextValueIndex];
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
		for (short index = 0; index < nextValueIndex; index++) {
			map.put(indexMap.get(index), countArray[index]);
		}
		map.put(null, entries);

		return map;
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

	public boolean isWithCounts() {
		return withCounts;
	}

	@Override
	public boolean put(CGATRingBuffer buffer, V value, boolean reverse) {
		long kmer = reverse ? CGAT.kMerToLongReverse(buffer) : CGAT.kMerToLongStraight(buffer);
		return putLong(kmer, value);
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
	 * 
	 * @return true if put succeeded, false if (probably) some value is already
	 *         stored under that kmer.
	 */
	protected boolean putLong(long kmer, V value) {
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
			if (nextValueIndex == Short.MAX_VALUE) {
				throw new IllegalStateException("Too many different values - only " + MAX_VALUES + " are possible.");
			}
			valueMap.put(value, nextValueIndex);
			indexMap.put(nextValueIndex, value);
			sindex = nextValueIndex;
			nextValueIndex++;
		}
		return sindex;
	}

	@Override
	public V get(CGATRingBuffer buffer, boolean reverse) {
		long kmer = reverse ? CGAT.kMerToLongReverse(buffer) : CGAT.kMerToLongStraight(buffer);
		return getLong(kmer, null);
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

	public byte getCount(long index) {
		if (largeCounts != null) {
			return BigArrays.get(largeCounts, index);
		} else {
			return counts[(int) index];
		}
	}

	public byte incCount(long index) {
		byte count = getCount(index);
		if (count == EXCLUDED_KMER_VIA_COUNT) {
			return count;
		}
		countsChanged = true;
		if (largeCounts != null) {
			BigArrays.incr(largeCounts, index);
			return ++count;
		} else {
			return ++counts[(int) index];
		}
	}

	public boolean update(CGATRingBuffer buffer, UpdateValueProvider<V> provider, boolean reverse) {
		if (!sorted) {
			throw new IllegalStateException("Updated only works when optimized.");
		}

		long kmer = reverse ? CGAT.kMerToLongReverse(buffer) : CGAT.kMerToLongStraight(buffer);

		long pos;
		short index;
		if (largeKmers != null) {
			pos = LongBigArrays.binarySearch(largeKmers, 0, entries, kmer);
			if (pos < 0) {
				return false;
			}
			index = BigArrays.get(largeValueIndexes, pos);
		} else {
			pos = Arrays.binarySearch(kmers, 0, (int) entries, kmer);
			if (pos < 0) {
				return false;
			}
			index = valueIndexes[(int) pos];
		}
		V oldValue = indexMap.get(index);
		V newValue = provider.getUpdateValue(oldValue);
		if (newValue == null) {
			throw new NullPointerException("null is not allowed as a value.");
		}
		if (newValue != oldValue && !newValue.equals(oldValue)) {
			if (largeKmers != null) {
				BigArrays.set(largeValueIndexes, pos, index);
			} else {
				valueIndexes[(int) pos] = index;
			}
		}
		return true;
	}

	public V getLong(long kmer, long[] posStore) {
		if (filter != null && !filter.containsLong(kmer)) {
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
		if (withCounts && getCount(index) == EXCLUDED_KMER_VIA_COUNT) {
			return null;
		}
		return indexMap.get(index);
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
					if (withCounts && countsChanged) {
						byte countA = BigArrays.get(largeCounts, a);
						byte countB = BigArrays.get(largeCounts, b);
						BigArrays.set(largeCounts, b, countA);
						BigArrays.set(largeCounts, a, countB);
					}
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
					if (withCounts && countsChanged) {
						byte countA = counts[(int) a];
						byte countB = counts[(int) b];
						counts[(int) b] = countA;
						counts[(int) a] = countB;
					}
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
			value = indexMap.get(sindex);
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
			if (!withCounts || getCount(i) != EXCLUDED_KMER_VIA_COUNT) {
				visitor.nextValue(this, kmer, sindex, i);
			}
		}
	}

	public interface KMerSortedArrayVisitor<V extends Serializable> {
		public void nextValue(KMerSortedArray<V> trie, long kmer, short index, long i);
	}

	public interface UpdateValueProvider<V extends Serializable> {
		public V getUpdateValue(V oldValue);
	}
}
