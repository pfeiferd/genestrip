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

import java.io.File;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;
import java.util.Iterator;

import org.metagene.genestrip.io.StreamProvider;

import it.unimi.dsi.fastutil.objects.Object2LongMap;

/**
 * A store mapping k-mers (each encoded as a {@code long}) to a value of type {@code V}, keyed
 * additionally by a compact per-store value index. Implementations are filled via {@link #putLong},
 * finalized once via {@link #optimize()}, and then serve {@link #getLong} lookups.
 *
 * @param <V> the value type mapped to each k-mer
 */
public interface KMerStore<V extends Serializable> extends Serializable {
	/**
	 * Reserves storage for the given number of expected k-mer entries. Must be called once before any
	 * {@link #putLong} insertions.
	 *
	 * @param size the number of expected k-mer entries to reserve storage for
	 */
	public void initSize(long size);

	/**
	 * Returns the k-mer length (k) of this store.
	 *
	 * @return the k-mer length
	 */
	public int getK();

	/**
	 * Returns the number of k-mer entries currently stored.
	 *
	 * @return the number of stored entries
	 */
	public long getEntries();

	/**
	 * Returns the reserved capacity of this store.
	 *
	 * @return the reserved capacity
	 */
	public long getSize();

	/**
	 * Finalizes the store for lookups (e.g. sorts the entries and rebuilds the pre-filter), enabling
	 * fast {@link #getLong} lookups and {@link #update}.
	 */
	public void optimize();

	/**
	 * Returns whether this store has been optimized for lookups.
	 *
	 * @return {@code true} if {@link #optimize()} has been called
	 */
	public boolean isOptimized();

	// --- Value <-> store-index mapping ---------------------------------------

	/**
	 * Returns the value-index capacity of this store.
	 *
	 * @return the maximum number of distinct values this store can hold (the value-index capacity).
	 *         Defined per implementation; see e.g. {@link KMerSortedArray#MAX_VALUES} and
	 *         {@link RadixKMerStore#maxValuesForRadix(int)}.
	 */
	public int getMaxValues();

	/**
	 * Returns the number of distinct values currently registered.
	 *
	 * @return the number of distinct values currently registered in this store.
	 */
	public int getNValues();

	/**
	 * Returns the value registered under the given store index.
	 *
	 * @param index the store index to look up
	 * @return the value registered under the given store index, or {@code null}.
	 */
	public V getValueForIndex(int index);

	/**
	 * Returns the store index of the given value.
	 *
	 * @param value the value to look up
	 * @return the store index of the given value, or {@code -1} if it is not registered.
	 */
	public int getIndexForValue(V value);

	/**
	 * Registers the value if necessary and returns its store index.
	 *
	 * @param value the value to register and/or look up
	 * @return the store index of the given value
	 */
	public int getAddValueIndex(V value);

	/**
	 * Returns an iterator over all registered values.
	 *
	 * @return an iterator over all values currently registered in this store.
	 */
	public Iterator<V> getValues();

	// --- Primary (long-keyed) access -----------------------------------------

	/**
	 * Inserts a value under the given k-mer.
	 *
	 * @param kmer  the k-mer (encoded as a {@code long}) to store the value under
	 * @param value the value to store
	 * @return {@code true} if the put succeeded, {@code false} if (probably) a value
	 *         is already stored under that k-mer.
	 */
	public boolean putLong(long kmer, V value);

	/**
	 * Looks up the value stored under the given k-mer.
	 *
	 * @param kmer     the k-mer (encoded as a {@code long}) to look up
	 * @param posStore optional single-element array; if non-{@code null}, the storage
	 *                 position of the k-mer is written to {@code posStore[0]}.
	 * @return the stored value, or {@code null} if the k-mer is not present.
	 */
	public V getLong(long kmer, long[] posStore);

	/**
	 * Returns whether the store is full.
	 *
	 * @return {@code true} if the store has reached its initialized capacity.
	 */
	public boolean isFull();

	/**
	 * Updates the value stored under the given k-mer using the given provider.
	 *
	 * @param kmer     the k-mer (encoded as a {@code long}) whose value is updated
	 * @param provider supplies the new value from the currently stored value
	 * @return {@code true} if the stored value was changed.
	 */
	public boolean update(long kmer, UpdateValueProvider<V> provider);

	/**
	 * Visits all stored entries, exposing the store index and storage position of each.
	 *
	 * @param visitor the visitor notified for each stored entry
	 */
	public void visit(IndexedKMerStoreVisitor<V> visitor);

	/**
	 * Stores the given value index for the entry at the given storage position. The position is one
	 * previously reported via {@link #getLong(long, long[])}'s {@code posStore} or
	 * {@link IndexedKMerStoreVisitor#nextValue(KMerStore, long, int, long)}.
	 *
	 * @param pos   the storage position of the entry
	 * @param index the value index to store
	 */
	public void setIndexAtPosition(long pos, int index);

	// --- Statistics -----------------------------------------------------------

	/**
	 * Fixes (caches) the number of k-mers per value so that it survives value updates.
	 */
	public void fix();

	/**
	 * Returns the number of stored k-mers per value.
	 *
	 * @return the (possibly cached, see {@link #fix()}) number of stored k-mers per value.
	 */
	public Object2LongMap<V> getFixedNKmersPerTaxid();

	/**
	 * Counts, by visiting every entry, how many stored k-mers map to each value. The {@code null} key
	 * maps to the total number of entries.
	 *
	 * @return a map from each value to its number of stored k-mers, with the {@code null} key
	 *         mapping to the total number of entries
	 */
	public Object2LongMap<V> getNKmersPerTaxid();

	// --- Value conversion -----------------------------------------------------

	/**
	 * Produces a new store of the same kind whose values are derived from this store's
	 * values via the given converter. The underlying k-mer data is shared where possible.
	 *
	 * @param <W>       the value type of the produced store
	 * @param converter converts each value of this store to a value of the new store
	 * @return a new store holding the converted values
	 */
	public <W extends Serializable> KMerStore<W> convertValues(ValueConverter<V, W> converter);

	/**
	 * Serializes the given store to the given file.
	 *
	 * @param <T>      the value type of the store
	 * @param store    the store to serialize
	 * @param trieFile the file to write the store to
	 * @throws IOException if writing the file fails
	 */
	public static <T extends Serializable> void save(KMerStore<T> store, File trieFile) throws IOException {
		try (ObjectOutputStream oOut = new ObjectOutputStream(StreamProvider.getOutputStreamForFile(trieFile))) {
			oOut.writeObject(store);			
		}
	}

	/**
	 * Deserializes a store previously written by {@link #save(KMerStore, File)}.
	 *
	 * @param <T>        the value type of the store
	 * @param filterFile the file to read the store from
	 * @return the deserialized store
	 * @throws IOException            if reading the file fails
	 * @throws ClassNotFoundException if the serialized store class cannot be found
	 */
	@SuppressWarnings("unchecked")
	public static <T extends Serializable> KMerStore<T> load(File filterFile)
			throws IOException, ClassNotFoundException {
		try (ObjectInputStream oOut = new ObjectInputStream(StreamProvider.getInputStreamForFile(filterFile))) {			
			return (KMerStore<T>) oOut.readObject();
		}
	}

	/**
	 * Visitor that exposes the store index of the value and the storage position of the entry,
	 * see {@link KMerStore#visit(IndexedKMerStoreVisitor)}.
	 *
	 * @param <V> the value type of the visited store
	 */
	public interface IndexedKMerStoreVisitor<V extends Serializable> {
		/**
		 * Called for each stored entry during a {@link KMerStore#visit(IndexedKMerStoreVisitor)}.
		 *
		 * @param store the store being visited
		 * @param kmer the k-mer (encoded as a {@code long}) of the entry
		 * @param index the store index of the entry's value
		 * @param pos the storage position of the entry
		 */
		public void nextValue(KMerStore<V> store, long kmer, int index, long pos);
	}

	/**
	 * Provides the new value for an {@link KMerStore#update(long, UpdateValueProvider)}.
	 *
	 * @param <V> the value type of the store being updated
	 */
	public interface UpdateValueProvider<V extends Serializable> {
		/**
		 * Returns the new value to store, derived from the currently stored value.
		 *
		 * @param oldValue the currently stored value
		 * @return the new value to store
		 */
		public V getUpdateValue(V oldValue);
	}

	/**
	 * Converts a value of type {@code V} into a value of type {@code W}, see
	 * {@link KMerStore#convertValues(ValueConverter)}.
	 *
	 * @param <V> the source value type
	 * @param <W> the target value type
	 */
	public interface ValueConverter<V, W extends Serializable> {
		/**
		 * Converts the given source value to a target value.
		 *
		 * @param value the source value to convert
		 * @return the converted target value
		 */
		public W convertValue(V value);
	}
}
