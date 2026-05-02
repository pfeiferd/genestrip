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
package org.metagene.genestrip.bloom;

import org.metagene.genestrip.util.LargeBitVector;

import java.util.Random;

public abstract class AbstractKMerBloomFilter implements KMerProbFilter {
	private static final long serialVersionUID = 1L;

	protected final double fpp;
	protected final Random random;

	protected long expectedInsertions;
	protected long bits;
	protected LargeBitVector bitVector;
	protected boolean largeBV;
	protected int hashes;
	protected long[] hashFactors;
	protected long entries;
	// For thread-safe put:
	protected long[] entriesAddends;
	protected final Object[] syncLocks;

	public AbstractKMerBloomFilter(double fpp) {
		if (fpp <= 0 || fpp >= 1) {
			throw new IllegalArgumentException("fpp must be a probability");
		}
		this.fpp = fpp;

		random = new Random(42);

		bitVector = new LargeBitVector(0);
		largeBV = bitVector.isLarge();
		entries = 0;
		bits = 0;
		syncLocks = new Object[256];
		for (int i = 0; i < 256; i++) {
			syncLocks[i] = new Object();
		}
	}

	@Override
	public void clear() {
		bitVector.clear();
		entries = 0;
	}

	@Override
	public long ensureExpectedSize(long expectedInsertions, boolean enforceLarge) {
		if (expectedInsertions < 0) {
			throw new IllegalArgumentException("expected insertions must be > 0");
		}
		this.expectedInsertions = expectedInsertions;

		bits = optimalNumOfBits(expectedInsertions, fpp);
		if (bitVector.ensureCapacity(bits, enforceLarge)) {
			hashes = optimalNumOfHashFunctions(expectedInsertions, bits);
			hashFactors = new long[hashes];
			for (int i = 0; i < hashFactors.length; i++) {
				hashFactors[i] = random.nextLong();
			}
		}
		largeBV = bitVector.isLarge();
		return bits;
	}

	public long getExpectedInsertions() {
		return expectedInsertions;
	}

	public int getHashes() {
		return hashes;
	}

	public double getFpp() {
		return fpp;
	}

	public long getBitSize() {
		return bitVector.getBitSize();
	}

	public boolean isLarge() {
		return bitVector.isLarge();
	}

	@Override
	public long getEntries() {
		if (entriesAddends != null) {
			long sum = entries;
			for (int i = 0; i < entriesAddends.length; i++) {
				sum += entriesAddends[i];
			}
			return sum;
		}
		else {
			return entries;
		}
	}

	protected int optimalNumOfHashFunctions(long n, long m) {
		return Math.max(1, (int) Math.round(((double) m) / n * Math.log(2)));
	}

	protected long optimalNumOfBits(long n, double p) {
		return (long) (-n * Math.log(p) / (Math.log(2) * Math.log(2)));
	}

	public void putLongThreadSafe(final long data) {
		for (int i = 0; i < hashes; i++) {
			long bitPos = reduce(hash(data, i));
			int syncIndex = ((int) ((byte) bitPos)) + 128;
			// Does this really bring a performance gain??
			// Locking is likely collision-free but must be done `hashes` times with computation of syncIndex...
			synchronized (syncLocks[syncIndex]) {
				bitVector.set(bitPos);
				if (i == 0) {
					if (entriesAddends == null) {
						entriesAddends = new long[256];
					}
					entriesAddends[syncIndex]++;
				}
			}
		}
	}

	@Override
	public void putLong(final long data) {
		entries++;
		for (int i = 0; i < hashes; i++) {
			bitVector.set(reduce(hash(data, i)));
		}
	}

	@Override
	public boolean containsLong(final long data) {
		for (int i = 0; i < hashes; i++) {
			if (!bitVector.get(reduce(hash(data, i)))) {
				return false;
			}
		}
		return true;
	}

	protected abstract long hash(final long data, final int i);

	protected long reduce(final long v) {
		return (v < 0 ? -v : v) % bits;
	}
}