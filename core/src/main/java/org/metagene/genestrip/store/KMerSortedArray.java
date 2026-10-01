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
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.concurrent.atomic.AtomicBoolean;

import org.metagene.genestrip.probfilter.BloomFilter;
import org.metagene.genestrip.probfilter.ProbFilter;
import org.metagene.genestrip.probfilter.MurmurBloomFilter;
import org.metagene.genestrip.probfilter.XORBloomFilter;
import org.metagene.genestrip.util.CGAT;
import org.metagene.genestrip.util.LargeBitVector;

import it.unimi.dsi.fastutil.BigArrays;
import it.unimi.dsi.fastutil.BigSwapper;
import it.unimi.dsi.fastutil.longs.LongBigArrays;
import it.unimi.dsi.fastutil.longs.LongComparator;

/**
 * A {@link KMerStore} that keeps its k-mers in a single {@code long} array sorted by
 * {@link #optimize()}, with a parallel {@code short} array holding each k-mer's value index (offset by
 * {@link Short#MIN_VALUE}); lookups binary-search the sorted array. It uses small ({@code int}-indexed)
 * arrays up to {@link #MAX_SMALL_CAPACITY} and switches to large fastutil
 * {@link it.unimi.dsi.fastutil.BigArrays} for bigger databases (or when large storage is enforced).
 * Common value-index, counter, pre-filter and statistics handling is inherited from
 * {@link AbstractKMerStore}.
 *
 * @param <V> the value type mapped to each k-mer
 */
public class KMerSortedArray<V extends Serializable> extends AbstractKMerStore<V> {
	// The value index is stored in a (large) short array, offset by Short.MIN_VALUE, so the number
	// of distinct values is capped by the short range.
	/** Maximum number of distinct values, capped by the short range of the value index. */
	public static final int MAX_VALUES = ((int) Short.MAX_VALUE) - Short.MIN_VALUE;

	private static final long serialVersionUID = 2L;

	/** Capacity threshold above which the store switches to big (paged) arrays. */
	public static final long MAX_SMALL_CAPACITY = Integer.MAX_VALUE - 8;

	/** The k-mers, sorted after optimization (small-array layout). */
	private long[] kmers;
	/** The per-k-mer value indexes (small-array layout). */
	private short[] valueIndexes;

	/** The k-mers, sorted after optimization (big-array layout). */
	private long[][] largeKmers;
	/** The per-k-mer value indexes (big-array layout). */
	private short[][] largeValueIndexes;

	/** Whether the big-array layout is enforced regardless of size. */
	private final boolean enforceLarge;

	// Just for optimizing synchronization during updates.
	/** The lock array used to serialize concurrent updates of the same entry. */
	private transient Object[] syncs;

	// Transient: marking state belongs to a run, never to a saved database.
	/** The visited marks, one bit per storage position, or {@code null} while marking was never on. */
	private transient LargeBitVector visited;
	/** Whether lookups through <em>this</em> instance mark the positions they hit; the hot-path flag. */
	private transient boolean markVisited;
	/**
	 * Whether any user of this store is marking. The marks live in the {@link #visited} vector, which
	 * a converted copy shares just as it shares the k-mer array, so two matchers marking at the same
	 * time would count each other's hits.
	 */
	private final AtomicBoolean markClaim = new AtomicBoolean(false);

	/**
	 * Creates a store with a fill-time pre-filter of the given {@code entryFpp} (an
	 * {@link XORBloomFilter} when {@code xor}, otherwise a {@link MurmurBloomFilter}) and a
	 * post-optimize filter targeting {@code optimizedFpp}, sized to hold {@code size} k-mers.
	 *
	 * @param k the k-mer length
	 * @param entryFpp the target false-positive probability of the fill-time pre-filter
	 * @param optimizedFpp the target false-positive probability after optimization
	 * @param initialValues the initial list of values
	 * @param enforceLarge whether to enforce the big-array layout
	 * @param xor whether to use an {@link XORBloomFilter} (else a {@link MurmurBloomFilter})
	 * @param size the number of expected k-mer entries to reserve storage for (must be {@code >= 0})
	 */
	/**
	 * Copies every entry of the given store into a new store of this class, i.e. from whatever layout
	 * the source uses into one sorted array searched by bisection. The copy is sized to the source's
	 * entry count and holds the same values, so the two differ in the data structure a lookup walks
	 * and in nothing else. This is what a comparison of the two needs: one database, two layouts, in
	 * one JVM.
	 * <p>
	 * The entries are written straight into the arrays, without the fill-time filter that
	 * {@link #putLong} consults. That filter exists to recognize a k-mer that is already stored, and a
	 * store being copied holds each of its k-mers once, so there is nothing to recognize. Leaving it
	 * out keeps every entry - {@code putLong} drops one whenever the filter answers a false positive -
	 * and saves its memory, which at the fill-time default of {@code 1e-11} is some 53 bits per k-mer,
	 * i.e. more than half of what the copied arrays themselves take.
	 *
	 * @param <V> the value type of both stores
	 * @param source the store to copy, which is left untouched and must be optimized
	 * @param optimizedFpp the target false-positive probability of the copy after optimization
	 * @return the copy, optimized and ready for lookups
	 * @throws IllegalArgumentException if the source holds more values than {@link #MAX_VALUES}
	 */
	public static <V extends Serializable> KMerSortedArray<V> copyOf(KMerStore<V> source, double optimizedFpp) {
		List<V> values = new ArrayList<V>();
		for (Iterator<V> i = source.getValues(); i.hasNext();) {
			values.add(i.next());
		}
		if (values.size() > MAX_VALUES) {
			throw new IllegalArgumentException("The source store holds " + values.size()
					+ " values, but this store caps them at " + MAX_VALUES + ".");
		}
		final KMerSortedArray<V> target = new KMerSortedArray<V>(source.getK(), values, false, null, optimizedFpp,
				source.getEntries());
		// The value indices of the two stores need not agree, so they are mapped on first use. The
		// map is as long as the source's index space, which is a megabyte at most.
		final int[] valueIndexMap = new int[source.getMaxValues()];
		Arrays.fill(valueIndexMap, -1);
		source.visit(new KMerStore.IndexedKMerStoreVisitor<V>() {
			@Override
			public void nextValue(KMerStore<V> store, long kmer, int index, long pos) {
				int mapped = valueIndexMap[index];
				if (mapped < 0) {
					mapped = target.getAddValueIndex(store.getValueForIndex(index));
					valueIndexMap[index] = mapped;
				}
				target.append(kmer, mapped);
			}
		});
		target.optimize();
		target.setUseFilter(source.isUseFilter() && target.getFilter() != null);
		return target;
	}

	/**
	 * Appends one entry, for a bulk copy from another store: the caller has already established that
	 * the k-mer is not stored yet, so no filter is consulted and no lock is taken. Single-threaded by
	 * contract, as {@link KMerStore#visit(IndexedKMerStoreVisitor)} is.
	 *
	 * @param kmer the k-mer to store
	 * @param valueIndex the index of its value in this store
	 */
	private void append(long kmer, int valueIndex) {
		if (entries == size) {
			throw new IllegalStateException("The copy was sized for " + size + " entries.");
		}
		long pos = entries++;
		if (largeKmers != null) {
			BigArrays.set(largeKmers, pos, kmer);
			BigArrays.set(largeValueIndexes, pos, (short) (valueIndex + Short.MIN_VALUE));
		} else {
			kmers[(int) pos] = kmer;
			valueIndexes[(int) pos] = (short) (valueIndex + Short.MIN_VALUE);
		}
		sorted = false;
	}

	public KMerSortedArray(int k, double entryFpp, double optimizedFpp, List<V> initialValues, boolean enforceLarge, boolean xor, long size) {
		this(k, initialValues, enforceLarge,
				xor ? new XORBloomFilter(entryFpp, size) : new BloomFilter(entryFpp, size), optimizedFpp, size);
	}

	/**
	 * Creates a store using the given fill-time pre-filter, sized to hold {@code size} k-mers. The
	 * backing k-mer/value arrays are allocated here and the pre-filter is expected to already be sized
	 * for {@code size} insertions, so the store is ready for {@link #putLong} insertions immediately
	 * after construction.
	 *
	 * @param k the k-mer length
	 * @param initialValues the initial list of values
	 * @param enforceLarge whether to enforce the big-array layout
	 * @param filter the fill-time pre-filter to use, already sized for {@code size} insertions
	 * @param optimizedFpp the target false-positive probability after optimization
	 * @param size the number of expected k-mer entries to reserve storage for (must be {@code >= 0})
	 */
	protected KMerSortedArray(int k, List<V> initialValues, boolean enforceLarge,
							  ProbFilter filter, double optimizedFpp, long size) {
		super(k, MAX_VALUES, initialValues, filter, optimizedFpp);
		this.enforceLarge = enforceLarge;
		if (size < 0) {
			throw new IllegalArgumentException("Expected insertions must be >= 0.");
		}
		this.size = size;
		if (size > MAX_SMALL_CAPACITY || enforceLarge) {
			largeKmers = BigArrays.ensureCapacity(BigArrays.wrap(new long[0]), size);
			largeValueIndexes = BigArrays.ensureCapacity(BigArrays.wrap(new short[0]), size);
		} else {
			kmers = new long[(int) size];
			valueIndexes = new short[(int) size];
		}
		initSyncs();
	}

	/**
	 * Value-converting copy constructor: shares the (immutable after optimize) k-mer and value-index
	 * arrays with {@code org} while remapping the values via {@code converter}.
	 *
	 * @param <W> the source value type
	 * @param org the store to copy k-mer and value-index arrays from
	 * @param converter the converter remapping source values to this store's value type
	 */
	public <W extends Serializable> KMerSortedArray(KMerSortedArray<W> org, KMerStore.ValueConverter<W, V> converter) {
		super(org, converter);
		kmers = org.kmers;
		valueIndexes = org.valueIndexes;
		largeKmers = org.largeKmers;
		largeValueIndexes = org.largeValueIndexes;
		visited = org.visited;
		enforceLarge = org.enforceLarge;
		initSyncs();
	}

	@Override
	public <W extends Serializable> KMerStore<W> convertValues(KMerStore.ValueConverter<V, W> converter) {
		return new KMerSortedArray<>(this, converter);
	}

	/**
	 * Indicates whether this store uses the big-array layout for large k-mer sets.
	 *
	 * @return whether this store uses the big-array layout
	 */
	public boolean isLarge() {
		return largeKmers != null;
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
				// Overfull store, leave quietly.
				return false;
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

	/**
	 * Decodes the k-mer starting at {@code start} in the nucleotide sequence and returns its stored
	 * value, see {@link #getLong(long, long[])}.
	 *
	 * @param nseq the nucleotide sequence
	 * @param start the start position of the k-mer within the sequence
	 * @param indexStore an optional array receiving the storage position (may be {@code null})
	 * @return the stored value, or {@code null} if the k-mer is not present
	 */
	public V get(byte[] nseq, int start, long[] indexStore) {
		return getLong(CGAT.kMerToLong(nseq, start, k, null), indexStore);
	}

	/**
	 * Restores the transient stripe-lock pool after deserialization.
	 *
	 * @param ois the stream to read from
	 * @throws ClassNotFoundException if a serialized class cannot be resolved
	 * @throws IOException            if reading fails
	 */
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

	/**
	 * Returns the stripe lock for the given position, mapping it onto the fixed {@link #syncs} pool with
	 * a cheap mask (the pool size is a power of two, so no modulo is needed): the same position always
	 * yields the same lock, different positions usually different ones, so updates only contend when
	 * they might touch the same entry. The key is a {@code long} so that a big-array position cannot
	 * overflow the mixing.
	 *
	 * @param key the storage position to derive the lock from
	 * @return the stripe lock to synchronize on
	 */
	private Object syncFor(long key) {
		return syncs[(int) (key & (syncs.length - 1))];
	}

	@Override
	public boolean update(long kmer, KMerStore.UpdateValueProvider<V> provider) {
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
		// We only must synchronize if two threads access the same kmer in the same position. Since we
		// cannot synchronize on the primitive 'pos', syncFor(pos) maps it onto a fixed stripe-lock pool:
		// the same pos always yields the same lock, different positions usually different ones. This
		// trick greatly decreases synchronization bottlenecks.
		synchronized (syncFor(pos)) {
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
				// Values are all pre-registered before the update phase, so this is a lock-free read.
				index = getRegisteredValueIndex(newValue);
				setIndexAtPosition(pos, index);
				return true;
			}
			return false;
		}
	}

	/**
	 * Stores the given value index for the entry at the given storage position.
	 *
	 * @param pos the storage position
	 * @param index the value index to store
	 */
	@Override
	public void setIndexAtPosition(long pos, int index) {
		if (largeKmers != null) {
			BigArrays.set(largeValueIndexes, pos, (short) (index + Short.MIN_VALUE));
		} else {
			valueIndexes[(int) pos] = (short) (index + Short.MIN_VALUE);
		}
		// This entry's value changed, so any cached per-value k-mer counts are stale.
		invalidateNKmersPerTaxid();
	}

	/**
	 * Returns the k-mer stored at the given storage position.
	 *
	 * @param pos the storage position
	 * @return the k-mer stored at the given storage position.
	 */
	public long getKMerAt(long pos) {
		if (largeKmers != null) {
			return BigArrays.get(largeKmers, pos);
		}
		return kmers[(int) pos];
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
		if (markVisited) {
			visited.set(pos);
		}
		return indexMap[index - Short.MIN_VALUE];
	}

	/**
	 * Returns the value index of the entry at the given storage position.
	 *
	 * @param pos the storage position
	 * @return the value index of the entry at the given storage position.
	 */
	public int indexAtPosition(long pos) {
		return (largeKmers != null ? BigArrays.get(largeValueIndexes, pos) : valueIndexes[(int) pos]) - Short.MIN_VALUE;
	}

	/**
	 * Turns the "visited" marking of lookups on or off. While it is on, every lookup that finds a
	 * k-mer sets the bit of its storage position, so the distinct k-mers a run has matched can be
	 * counted afterwards with {@link #countVisitedPerValueIndex(long[])}.
	 * <p>
	 * Unlike {@link RadixKMerStore}, which marks the entry it has just loaded, the marks live in a
	 * bit vector of their own: one more random access per matched k-mer, and about 1.6 % of the k-mer
	 * array in memory. Only one user can mark at a time, since a converted copy shares the vector.
	 *
	 * @param markVisited whether lookups should mark the positions they hit
	 * @return whether marking could be enabled; {@code false} means someone else is marking
	 */
	@Override
	public boolean setMarkVisited(boolean markVisited) {
		if (markVisited) {
			// Re-entrant for the holder: a matcher claims once per run and may start several runs.
			if (!this.markVisited && !markClaim.compareAndSet(false, true)) {
				return false;
			}
			synchronized (this) {
				if (visited == null) {
					visited = new LargeBitVector(entries);
				}
			}
		} else if (this.markVisited) {
			markClaim.set(false);
		}
		this.markVisited = markVisited;
		return true;
	}

	@Override
	public boolean isMarkVisited() {
		return markVisited;
	}

	@Override
	public void clearVisitedMarks() {
		if (visited != null) {
			visited.clear();
		}
	}

	@Override
	public void countVisitedPerValueIndex(long[] counts) {
		if (visited == null) {
			return;
		}
		for (long i = 0; i < entries; i++) {
			if (visited.get(i)) {
				counts[indexAtPosition(i)]++;
			}
		}
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
		// Rework the bloom filter from the fill-time one into the optimized one.
		filter = createOptimizedFilter();
		if (filter != null) {
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
	public void visit(KMerStore.IndexedKMerStoreVisitor<V> visitor) {
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
}
