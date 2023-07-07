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
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.metagene.genestrip.bloom.MurmurCGATBloomFilter;
import org.metagene.genestrip.util.CGAT;
import org.metagene.genestrip.util.CGATRingBuffer;

import it.unimi.dsi.fastutil.BigArrays;
import it.unimi.dsi.fastutil.BigSwapper;
import it.unimi.dsi.fastutil.longs.LongBigArrays;
import it.unimi.dsi.fastutil.longs.LongComparator;
import it.unimi.dsi.fastutil.objects.Object2ShortMap;
import it.unimi.dsi.fastutil.objects.Object2ShortOpenHashMap;
import it.unimi.dsi.fastutil.shorts.Short2ObjectMap;
import it.unimi.dsi.fastutil.shorts.Short2ObjectOpenHashMap;

public class KMerSortedArray<V extends Serializable> implements KMerStore<V> {
	private static final long serialVersionUID = 1L;

	public static long MAX_SMALL_CAPACITY = Integer.MAX_VALUE - 8;

	private long[] kmers;
	private short[] valueIndexes;

	private long[][] largeKmers;
	private short[][] largeValueIndexes;

	private final int k;
	private final Object2ShortMap<V> valueMap;
	private final Short2ObjectMap<V> indexMap;

	protected long size;

	private long entries;
	private boolean sorted;
	private boolean large;
	private boolean initSize;
	private short nextValueIndex;
	private MurmurCGATBloomFilter filter;

	public KMerSortedArray(int k, double fpp, List<V> initialValues) {
		this(k, fpp, initialValues, false);
	}

	public KMerSortedArray(int k, double fpp, List<V> initialValues, boolean enforceLarge) {
		this.k = k;
		int s = initialValues == null ? 0 : initialValues.size();
		indexMap = new Short2ObjectOpenHashMap<V>(s);
		valueMap = new Object2ShortOpenHashMap<V>(s);
		nextValueIndex = 1;
		if (initialValues != null) {
			for (V v : initialValues) {
				if (nextValueIndex == 0) {
					throw new IllegalStateException("Too many different values - only 65534 are possible.");
				}
				valueMap.put(v, nextValueIndex);
				indexMap.put(nextValueIndex, v);
				nextValueIndex++;
			}
		}
		filter = new MurmurCGATBloomFilter(k, fpp);
		large = enforceLarge;
	}

	@Override
	public void initSize(long expectedInsertions) {
		if (initSize) {
			throw new IllegalStateException("Cant initlialize size twice.");
		}
		initSize = true;
		clearAndEnsureCapacity(expectedInsertions);
	}

	public void clearAndEnsureCapacity(long expectedInsertions) {
		if (expectedInsertions <= 0) {
			throw new IllegalArgumentException("Expected insertions must be > 0.");
		}
		filter.clearAndEnsureCapacity(expectedInsertions);

		if (size >= expectedInsertions) {
			clearArray();
		} else {
			size = expectedInsertions;
			initBitArray();
		}
	}

	protected void initBitArray() {
		if (size > MAX_SMALL_CAPACITY || large == true) {
			large = true;
			kmers = null;
			valueIndexes = null;
			if (largeKmers == null) {
				largeKmers = BigArrays.ensureCapacity(BigArrays.wrap(new long[0]), size);
				largeValueIndexes = BigArrays.ensureCapacity(BigArrays.wrap(new short[0]), size);
			} else {
				largeKmers = BigArrays.ensureCapacity(largeKmers, size);
				largeValueIndexes = BigArrays.ensureCapacity(largeValueIndexes, size);
			}
		} else {
			kmers = new long[(int) size];
			valueIndexes = new short[(int) size];
		}
	}

	protected void clearArray() {
		if (large) {
			BigArrays.fill(largeKmers, 0);
			BigArrays.fill(largeValueIndexes, (short) 0);
		} else {
			Arrays.fill(kmers, 0);
			Arrays.fill(valueIndexes, (short) 0);
		}
	}

	@Override
	public Map<V, Long> getNKmersPerTaxid() {
		long[] countArray = new long[nextValueIndex - 1];
		if (large) {
			for (long i = 0; i < entries; i++) {
				short index = BigArrays.get(largeValueIndexes, i);
				countArray[index - 1]++;
			}
		} else {
			for (int i = 0; i < entries; i++) {
				short index = valueIndexes[i];
				countArray[index - 1]++;
			}
		}

		Map<V, Long> map = new HashMap<V, Long>();
		for (short index = 1; index < nextValueIndex; index++) {
			map.put(indexMap.get(index), countArray[index - 1]);
		}
		map.put(null, entries);

		return map;
	}

	@Override
	public int getMaxValues() {
		return 65534;
	}

	@Override
	public int getK() {
		return k;
	}

	@Override
	public long getEntries() {
		return entries;
	}

	public long getSize() {
		return size;
	}

	public boolean isLarge() {
		return large;
	}

	@Override
	public boolean put(CGATRingBuffer buffer, V value, boolean reverse) {
		long kmer = reverse ? CGAT.kmerToLongReverse(buffer) : CGAT.kmerToLongStraight(buffer);
		return putLong(kmer, value);
	}

	@Override
	public boolean put(byte[] nseq, int start, V value, boolean reverse) {
		long kmer = reverse ? CGAT.kMerToLongReverse(nseq, start, k, null)
				: CGAT.kMerToLongStraight(nseq, start, k, null);
		return putLong(kmer, value);
	}

	protected boolean putLong(long kmer, V value) {
		if (value == null) {
			throw new NullPointerException("null is not allowed as a value.");
		}
		sorted = false;
		if (entries == size) {
			throw new IllegalStateException("Capacity exceeded.");
		}
		if (filter.containsLong(kmer)) {
			V v = getLong(kmer, null);
			if (v != null) {
				return v.equals(value);
			}
		}
		short sindex = valueMap.getShort(value);
		if (sindex == valueMap.defaultReturnValue()) {
			if (nextValueIndex == 0) {
				throw new IllegalStateException("Too many different values - only 65534 are possible.");
			}
			valueMap.put(value, nextValueIndex);
			indexMap.put(nextValueIndex, value);
			sindex = nextValueIndex;
			nextValueIndex++;
		}
		if (large) {
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

	@Override
	public V get(CGATRingBuffer buffer, boolean reverse) {
		long kmer = reverse ? CGAT.kmerToLongReverse(buffer) : CGAT.kmerToLongStraight(buffer);
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

	protected V getLong(long kmer, long[] indexStore) {
		if (filter != null && !filter.containsLong(kmer)) {
			return null;
		}
		short index;
		if (sorted) {
			if (large) {
				long pos = LongBigArrays.binarySearch(largeKmers, 0, entries, kmer);
				if (pos < 0) {
					return null;
				}
				index = BigArrays.get(largeValueIndexes, pos);
			} else {
				int pos = Arrays.binarySearch(kmers, 0, (int) entries, kmer);
				if (pos < 0) {
					return null;
				}
				index = valueIndexes[pos];
			}
		} else {
			if (large) {
				long pos = -1;
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
				int pos = -1;
				for (int i = 0; i < entries; i++) {
					if (kmers[i] == kmer) {
						pos = i;
						break;
					}
				}
				if (pos < 0) {
					return null;
				}
				index = valueIndexes[pos];
			}
		}
		if (indexStore != null) {
			indexStore[0] = index;
		}
		return indexMap.get(index);
	}

	@Override
	public void optimize() {
		if (sorted) {
			return;
		}
		if (large) {
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
			if (large) {
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
}
