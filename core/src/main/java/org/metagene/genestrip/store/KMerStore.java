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

public interface KMerStore<V extends Serializable> extends Serializable {
	public void initSize(long size);

	public int getK();

	public long getEntries();

	public long getSize();

	public void optimize();

	public boolean isOptimized();

	// --- Value <-> store-index mapping ---------------------------------------

	/**
	 * @return the number of distinct values currently registered in this store.
	 */
	public int getNValues();

	/**
	 * @return the value registered under the given store index, or {@code null}.
	 */
	public V getValueForIndex(int index);

	/**
	 * @return the store index of the given value, or {@code -1} if it is not registered.
	 */
	public int getIndexForValue(V value);

	/**
	 * Registers the value if necessary and returns its store index.
	 */
	public int getAddValueIndex(V value);

	/**
	 * @return an iterator over all values currently registered in this store.
	 */
	public Iterator<V> getValues();

	// --- Primary (long-keyed) access -----------------------------------------

	/**
	 * Inserts a value under the given k-mer.
	 *
	 * @return {@code true} if the put succeeded, {@code false} if (probably) a value
	 *         is already stored under that k-mer.
	 */
	public boolean putLong(long kmer, V value);

	/**
	 * Looks up the value stored under the given k-mer.
	 *
	 * @param posStore optional single-element array; if non-{@code null}, the storage
	 *                 position of the k-mer is written to {@code posStore[0]}.
	 * @return the stored value, or {@code null} if the k-mer is not present.
	 */
	public V getLong(long kmer, long[] posStore);

	/**
	 * @return {@code true} if the store has reached its initialized capacity.
	 */
	public boolean isFull();

	/**
	 * Updates the value stored under the given k-mer using the given provider.
	 *
	 * @return {@code true} if the stored value was changed.
	 */
	public boolean update(long kmer, UpdateValueProvider<V> provider);

	/**
	 * @return the number of k-mers whose value was moved (changed) via {@link #update}.
	 *         The default implementation returns {@code 0} for stores that do not track this.
	 */
	default long getKMersMoved() {
		return 0;
	}

	/**
	 * Visits all stored entries, exposing the store index and storage position of each.
	 */
	public void visit(IndexedKMerStoreVisitor<V> visitor);

	// --- Statistics -----------------------------------------------------------

	/**
	 * Fixes (caches) the number of k-mers per value so that it survives value updates.
	 */
	public void fix();

	/**
	 * @return the (possibly cached, see {@link #fix()}) number of stored k-mers per value.
	 */
	public Object2LongMap<V> getFixedNKmersPerTaxid();

	// --- Value conversion -----------------------------------------------------

	/**
	 * Produces a new store of the same kind whose values are derived from this store's
	 * values via the given converter. The underlying k-mer data is shared where possible.
	 */
	public <W extends Serializable> KMerStore<W> convertValues(ValueConverter<V, W> converter);

	public static <T extends Serializable> void save(KMerStore<T> store, File trieFile) throws IOException {
		try (ObjectOutputStream oOut = new ObjectOutputStream(StreamProvider.getOutputStreamForFile(trieFile))) {
			oOut.writeObject(store);			
		}
	}

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
	 */
	public interface IndexedKMerStoreVisitor<V extends Serializable> {
		public void nextValue(KMerStore<V> store, long kmer, int index, long pos);
	}

	/**
	 * Provides the new value for an {@link KMerStore#update(long, UpdateValueProvider)}.
	 */
	public interface UpdateValueProvider<V extends Serializable> {
		public V getUpdateValue(V oldValue);
	}

	/**
	 * Converts a value of type {@code V} into a value of type {@code W}, see
	 * {@link KMerStore#convertValues(ValueConverter)}.
	 */
	public interface ValueConverter<V, W extends Serializable> {
		public W convertValue(V value);
	}
}
