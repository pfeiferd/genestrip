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

/**
 * Base class for {@link KMerProbFilter} implementations backed by a {@link LargeBitVector}. It sizes
 * the bit vector and the number of hash functions from the expected insertions and target
 * false-positive probability, and derives each of the {@code hashes} bit indices from a
 * subclass-supplied {@link #hash} function. Subclasses only provide the hash.
 */
public abstract class AbstractKMerBloomFilter implements KMerProbFilter {
	private static final long serialVersionUID = 2L;

	/** The target false-positive probability. */
	protected final double fpp;
	/** Source of randomness used to derive the hash factors. */
	protected final Random random;

	/** The expected number of insertions the filter is sized for. */
	protected long expectedInsertions;
	/** The number of bits in the backing bit vector. */
	protected long bits;
	/** The backing bit vector storing the filter's bits. */
	protected LargeBitVector bitVector;
	/** The number of hash functions applied per k-mer. */
	protected int hashes;
	/** The random factors used to derive the individual hashes. */
	protected long[] hashFactors;

	/**
	 * Creates a filter with the given target false-positive probability (which must lie strictly
	 * between 0 and 1), sized for {@code expectedInsertions} k-mers. The bit vector and the number of
	 * hash functions are derived from the sizing here, so the filter is ready for insertions
	 * immediately; its size is fixed for the lifetime of the filter (the backing bit vector is always
	 * the bucketed {@link LargeBitVector}, which itself grows only if needed, but this filter never
	 * re-derives its {@link #bits}/{@link #hashes} sizing).
	 *
	 * @param fpp the target false-positive probability, strictly between 0 and 1
	 * @param expectedInsertions the expected number of k-mers to be inserted (must be {@code >= 0})
	 */
	public AbstractKMerBloomFilter(double fpp, long expectedInsertions) {
		if (fpp <= 0 || fpp >= 1) {
			throw new IllegalArgumentException("fpp must be a probability");
		}
		if (expectedInsertions < 0) {
			throw new IllegalArgumentException("expected insertions must be >= 0");
		}
		this.fpp = fpp;
		random = new Random(42);

		this.expectedInsertions = expectedInsertions;
		bits = optimalNumOfBits(expectedInsertions, fpp);
		bitVector = new LargeBitVector(bits);
		hashes = optimalNumOfHashFunctions(expectedInsertions, bits);
		hashFactors = new long[hashes];
		for (int i = 0; i < hashFactors.length; i++) {
			hashFactors[i] = random.nextLong();
		}
	}

	@Override
	public void clear() {
		bitVector.clear();
	}

	/**
	 * Returns the expected number of insertions the filter is currently sized for.
	 *
	 * @return the expected number of insertions the filter is currently sized for.
	 */
	public long getExpectedInsertions() {
		return expectedInsertions;
	}

	/**
	 * Returns the target false-positive probability.
	 *
	 * @return the target false-positive probability.
	 */
	public double getFpp() {
		return fpp;
	}

	/**
	 * Returns the number of bits in the backing bit vector.
	 *
	 * @return the number of bits in the backing bit vector.
	 */
	public long getBitSize() {
		return bitVector.getBitSize();
	}

	/**
	 * Computes the optimal number of hash functions for the given sizing.
	 *
	 * @param n the expected number of insertions
	 * @param m the number of bits in the backing bit vector
	 * @return the optimal number of hash functions for {@code n} expected insertions into {@code m} bits.
	 */
	protected int optimalNumOfHashFunctions(long n, long m) {
		return Math.max(1, (int) Math.round(((double) m) / n * Math.log(2)));
	}

	/**
	 * Computes the optimal number of bits for the given sizing and false-positive probability.
	 *
	 * @param n the expected number of insertions
	 * @param p the target false-positive probability
	 * @return the optimal number of bits for {@code n} expected insertions at false-positive probability {@code p}.
	 */
	protected long optimalNumOfBits(long n, double p) {
		return Math.max(1L, (long) (-n * Math.log(p) / (Math.log(2) * Math.log(2))));
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

	/**
	 * Adds the given k-mer to the filter and reports whether it was newly added, computing each hash
	 * only once and setting the corresponding bits atomically. This combines the effect of a {@link
	 * #containsLong(long)} check followed by {@link #putLong(long)} into a single hashing pass, and is
	 * safe to call concurrently from multiple threads (each bit is set with an atomic OR, so no
	 * concurrent update is lost and false negatives cannot occur). Because inserting an already
	 * present k-mer only re-sets bits that are already set, the resulting filter state is identical to
	 * a plain {@code if (!containsLong(data)) putLong(data)} sequence.
	 * <p>
	 * The returned "newly added" flag is exact under single-threaded use. Under concurrent use two
	 * threads inserting the same absent k-mer may both observe it as new; this never affects the
	 * filter's membership answers.
	 *
	 * @param data the k-mer, encoded as a {@code long}, to add
	 * @return {@code true} if the k-mer was not already present (at least one of its bits was newly
	 *         set), {@code false} if it was already present
	 */
	@Override
	public boolean putLong(final long data) {
		boolean added = false;
		for (int i = 0; i < hashes; i++) {
			// Every bit must be set even once we know the element is new, so do not short-circuit.
			if (bitVector.set(reduce(hash(data, i)))) {
				added = true;
			}
		}
		return added;
	}

	/**
	 * Computes the {@code i}-th hash of the given k-mer. The default is a Lemire-style bit-mixing hash
	 * (moved here from the former {@code LemireOptBloomFilter}); subclasses may override it with a
	 * cheaper or otherwise preferable hash (see {@link XORKMerBloomFilter}, {@link MurmurKMerBloomFilter}).
	 *
	 * @param data the k-mer, encoded as a {@code long}, to hash
	 * @param i    the index of the hash function to apply
	 * @return the {@code i}-th hash of the given k-mer.
	 */
	protected long hash(final long data, final int i) {
		long x = data + hashFactors[i];
		x = (x ^ (x >>> 33)) * 0xff51afd7ed558ccdL;
		x = (x ^ (x >>> 33)) * 0xc4ceb9fe1a85ec53L;
		x = x ^ (x >>> 33);
		return x;
	}

	/**
	 * Maps a hash value onto a valid bit index of the backing bit vector.
	 *
	 * @param v the hash value to reduce
	 * @return the bit index in {@code [0, bits)} for the given hash value.
	 */
	protected long reduce(final long v) {
		return Math.abs(v % bits);
	}
}
