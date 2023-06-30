package org.metagene.genestrip.trie;

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
import it.unimi.dsi.fastutil.shorts.Short2ObjectLinkedOpenHashMap;
import it.unimi.dsi.fastutil.shorts.Short2ObjectMap;

public class KMerSortedArray<V extends Serializable> implements KMerStore<V> {
	public static long MAX_SMALL_CAPACITY = Integer.MAX_VALUE - 8;

	private long[] kmers;
	private short[] valueIndexes;

	private long[][] largeKmers;
	private short[][] largeValueIndexes;

	private final int k;
	private final Map<V, Short> valueMap;
	private final Short2ObjectMap<V> indexMap;

	protected long size;

	private long entries;
	private boolean sorted;
	private boolean large;
	private boolean initSize;
	private short nextValueIndex;
	private MurmurCGATBloomFilter filter;

	public KMerSortedArray(int k, double fpp, List<V> initialValues) {
		this.k = k;
		indexMap = new Short2ObjectLinkedOpenHashMap<V>(initialValues.size());
		valueMap = new HashMap<V, Short>(initialValues.size());
		for (V v : initialValues) {
			valueMap.put(v, nextValueIndex);
			indexMap.put(nextValueIndex, v);
			nextValueIndex++;
		}
		filter = new MurmurCGATBloomFilter(k, fpp);
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
	public int getLen() {
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
		long kmer = reverse ? CGAT.kmerToLongReverse(nseq, start, k, null)
				: CGAT.kmerToLongStraight(nseq, start, k, null);
		return putLong(kmer, value);
	}

	protected boolean putLong(long kmer, V value) {
		if (sorted) {
			throw new IllegalStateException("Cannot insert after optimize.");
		}
		if (entries == size) {
			throw new IllegalStateException("Capacity exceeded.");
		}
		if (filter.containsLong(kmer)) {
			return false;
		}
		Short index = valueMap.get(value);
		short sindex;
		if (index == null) {
			valueMap.put(value, nextValueIndex);
			indexMap.put(nextValueIndex, value);
			sindex = nextValueIndex;
			nextValueIndex++;
		} else {
			sindex = index;
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
		return getLong(kmer);
	}

	@Override
	public V get(byte[] nseq, int start, boolean reverse) {
		long kmer = reverse ? CGAT.kmerToLongReverse(nseq, start, k, null)
				: CGAT.kmerToLongStraight(nseq, start, k, null);
		return getLong(kmer);
	}

	protected V getLong(long kmer) {
		if (filter != null && !filter.containsLong(kmer)) {
			return null;
		}
		if (sorted) {
			short index;
			if (large) {
				long pos = LongBigArrays.binarySearch(largeKmers, kmer);
				if (pos < 0) {
					return null;
				}
				index = BigArrays.get(largeValueIndexes, pos);
			} else {
				int pos = Arrays.binarySearch(kmers, kmer);
				if (pos < 0) {
					return null;
				}
				index =valueIndexes[pos];
			}
			return indexMap.get(index);
		} else {
			throw new IllegalStateException("Get only works when optimized / sorted.");
		}
	}

	@Override
	public void optimize() {
		BigArrays.quickSort(0, entries, new LongComparator() {
			@Override
			public int compare(long k1, long k2) {
				return k1 < k2 ? -1 : k1 > k2 ? 1 : 0;
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
			}
			else {
				kmer = kmers[(int) i];
				sindex = valueIndexes[(int) i];
			}
			value = indexMap.get(sindex);
			visitor.nextValue(this, kmer, value);
		}
	}
}
