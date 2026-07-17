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

import java.util.HashMap;
import java.util.Map;

import org.metagene.genestrip.store.KMerStore.IndexedKMerStoreVisitor;
import org.metagene.genestrip.util.LargeBitVector;
import org.metagene.genestrip.util.LargeShortVector;

import it.unimi.dsi.fastutil.objects.Object2LongLinkedOpenHashMap;
import it.unimi.dsi.fastutil.objects.Object2LongMap;

/**
 * A {@link KMerUniqueCounter} that marks each matched k-mer in a {@link LargeBitVector} indexed by the
 * k-mer's storage position in the backing {@link KMerStore}, so every distinct k-mer is counted at
 * most once. When constructed with counts enabled, a parallel {@link LargeShortVector} also records
 * how often each k-mer was matched.
 */
public class KMerUniqueCounterBits implements KMerUniqueCounter {
	private static final int LOCKS = 512; // Must be a value in 2^n, n = 1,2,3,...
	private static final int LOCKS_MASK = LOCKS - 1;

	private final KMerStore<String> store;
	private final LargeBitVector bitVector;
	private final LargeShortVector countsVector;
	private final Object[] locks;

	/**
	 * Creates a counter over the given store. If {@code withCounts} is set, per-k-mer match counts are
	 * tracked in addition to the presence bits.
	 *
	 * @param store      the backing k-mer store
	 * @param withCounts if {@code true}, per-k-mer match counts are tracked as well
	 */
	public KMerUniqueCounterBits(KMerStore<String> store, boolean withCounts) {
		this.store = store;
		bitVector = new LargeBitVector(store.getEntries());
		countsVector = withCounts ? new LargeShortVector(store.getEntries()) : null;
		locks = new Object[LOCKS];
		for (int i = 0; i < locks.length; i++) {
			locks[i] = new Object();
		}
	}
	
	/**
	 * Returns whether per-k-mer match counts are being tracked.
	 *
	 * @return {@code true} if counting is enabled
	 */
	public boolean isWithCounts() {
		return countsVector != null;
	}

	@Override
	public void clear() {
		bitVector.clear();
		if (countsVector != null) {
			countsVector.clear();
		}
	}

	@Override
	public final synchronized void put(final long kmer, final String taxid, final long index) {
		bitVector.set(index);
		if (countsVector != null) {
			countsVector.inc(index);
		}
	}

	/**
	 * Marks the k-mer at the given storage position as matched (and increments its match count when
	 * counting is enabled), synchronizing per bit-array word to allow concurrent callers.
	 *
	 * @param index the storage position of the k-mer within the backing store
	 */
	public final void putInlined(final long index) {
		// LOCKS is a power of two, so the stripe is selected with a cheap AND (no modulo). Keying on the
		// word index (index >>> 6) makes concurrent writers to the same word share a lock.
		synchronized (locks[(int) ((index >>> 6) & LOCKS_MASK)]) {
			bitVector.set(index);
			if (countsVector != null) {
				countsVector.inc(index);
			}
		}
	}

	@Override
	public Object2LongMap<String> getUniqueKmerCounts() {
		Object2LongMap<String> res = new Object2LongLinkedOpenHashMap<String>();
		long[] valueCounter = new long[store.getNValues()];
		store.visit(new IndexedKMerStoreVisitor<String>() {
			@Override
			public void nextValue(KMerStore<String> trie, long kmer, int index, long i) {
				if (bitVector.get(i)) {
					valueCounter[index]++;
				}
			}
		});
		for (int i = 0; i < valueCounter.length; i++) {
			String taxid = store.getValueForIndex(i);
			res.put(taxid, valueCounter[i]);
		}

		return res;
	}

	/**
	 * Determines, per taxid, the {@code counts} highest per-k-mer match counts among that taxid's
	 * matched k-mers. The returned map also holds the overall top counts under the {@code null} key.
	 * Requires that this counter was created with counting enabled.
	 *
	 * @param counts the number of highest match counts to keep per taxid
	 * @return a map from taxid to its highest match counts, with the overall top counts under the {@code null} key
	 */
	public Map<String, short[]> getMaxCountsCounts(int counts) {
		Map<String, short[]> res = new HashMap<String, short[]>();
		short[] totalMaxCounts = new short[counts];
		res.put(null, totalMaxCounts);
		
		
		store.visit(new IndexedKMerStoreVisitor<String>() {
			@Override
			public void nextValue(KMerStore<String> trie, long kmer, int index, long i) {
				if (bitVector.get(i)) {
					String taxid = store.getValueForIndex(index);
					if (taxid != null) {						
						short[] target = res.get(taxid); 
						if (target == null) {
							target = new short[counts];
							res.put(taxid, target);
						}								
						short count = countsVector.get(i);
						updateMaxCounts(count, target);
						updateMaxCounts(count, totalMaxCounts);
					}
				}
			}
		});

		return res;
	}
	
	private void updateMaxCounts(short count, short[] target) {
		for (int j = 0; j < target.length; j++) {
			if (count > target[j]) {
				for (int k = target.length - 1; k > j; k--) {
					target[k] = target[k - 1];
				}
				target[j] = count;
				return;
			}
		}		
	}

	@Override
	public int getUniqueKmerCount(String taxid) {
		final int sindex = store.getIndexForValue(taxid);
		if (sindex < 0) {
			return 0;
		}
		int[] count = new int[1];
		store.visit(new IndexedKMerStoreVisitor<String>() {
			@Override
			public void nextValue(KMerStore<String> trie, long kmer, int index, long i) {
				if (index == sindex && bitVector.get(i)) {
					count[0]++;
				}
			}
		});

		return count[0];
	}

	/**
	 * Fills {@code target} (in descending order) with the highest per-k-mer match counts among the
	 * given taxid's matched k-mers.
	 *
	 * @param taxid  the taxid whose matched k-mers are inspected
	 * @param target the array to fill with the highest match counts in descending order
	 */
	public void getMaxCounts(String taxid, short[] target) {
		final int sindex = store.getIndexForValue(taxid);
		if (sindex < 0) {
			return;
		}
		store.visit(new IndexedKMerStoreVisitor<String>() {
			@Override
			public void nextValue(KMerStore<String> trie, long kmer, int index, long i) {
				if (index == sindex && bitVector.get(i)) {
					updateMaxCounts(countsVector.get(i), target);
				}
			}
		});
	}
}