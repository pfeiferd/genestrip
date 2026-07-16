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

import it.unimi.dsi.fastutil.BigArrays;
import org.metagene.genestrip.util.LargeBitVector;

import java.io.IOException;
import java.io.ObjectOutputStream;
import java.util.Random;

/**
 * Base class for {@link KMerProbFilter} implementations backed by a {@link LargeBitVector}. It sizes
 * the bit vector and the number of hash functions from the expected insertions and target
 * false-positive probability, and derives each of the {@code hashes} bit indices from a
 * subclass-supplied {@link #hash} function. Subclasses only provide the hash. The concurrent
 * entry-counting used by {@link #putLongIfAbsent(long)} is inherited from {@link
 * AbstractCountingKMerProbFilter}; the plain {@link #entries} counter is kept here to preserve the
 * serialized form.
 */
public abstract class AbstractKMerBloomFilter extends AbstractCountingKMerProbFilter {
	private static final long serialVersionUID = 1L;

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
	/** Whether the backing bit vector uses large storage. */
	protected boolean largeBV;
	/** The number of hash functions applied per k-mer. */
	protected int hashes;
	/** The random factors used to derive the individual hashes. */
	protected long[] hashFactors;
	/** The number of k-mers added via the plain {@link #putLong(long)} path. */
	protected long entries;

	/**
	 * Creates an empty filter with the given target false-positive probability (which must lie strictly
	 * between 0 and 1).
	 *
	 * @param fpp the target false-positive probability, strictly between 0 and 1
	 */
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
	}

	@Override
	public void clear() {
		bitVector.clear();
		entries = 0;
		if (concurrentEntries != null) {
			concurrentEntries.reset();
		}
	}

	@Override
	public long ensureExpectedSize(long expectedInsertions, boolean enforceLarge) {
		if (expectedInsertions < 0) {
			throw new IllegalArgumentException("expected insertions must be > 0");
		}
		ensureConcurrentEntries();
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

	/**
	 * Returns the expected number of insertions the filter is currently sized for.
	 *
	 * @return the expected number of insertions the filter is currently sized for.
	 */
	public long getExpectedInsertions() {
		return expectedInsertions;
	}

	/**
	 * Returns the number of hash functions applied per k-mer.
	 *
	 * @return the number of hash functions applied per k-mer.
	 */
	public int getHashes() {
		return hashes;
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
	 * Returns whether the backing bit vector uses large storage.
	 *
	 * @return whether the backing bit vector uses large storage.
	 */
	public boolean isLarge() {
		return bitVector.isLarge();
	}

	@Override
	public long getEntries() {
		return entries + concurrentEntryCount();
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
	public void putLong(final long data) {
		entries++;
		for (int i = 0; i < hashes; i++) {
            final long index = reduce(hash(data, i));
            if (bitVector.largeBits != null) {
                long arrayIndex = index >>> 6;
                BigArrays.set(bitVector.largeBits, arrayIndex, BigArrays.get(bitVector.largeBits, arrayIndex) | (1L << (index & 0b111111)));
            } else {
                // Using optimization instead of '%', see:
                // http://lemire.me/blog/2016/06/27/a-fast-alternative-to-the-modulo-reduction/
                // and Line 34 in https://github.com/FastFilter/fastfilter_java/blob/master/fastfilter/src/main/java/org/fastfilter/utils/Hash.java
                // Not sure whether it would also work for long - probably not.
                //int arrayIndex = (int) ((((index >>> 6) & 0xffffffffL) * (size & 0xffffffffL)) >>> 32);
                // Original code:
                int arrayIndex = (int) (index >>> 6);
                bitVector.bits[arrayIndex] = bitVector.bits[arrayIndex] | (1L << (index & 0b111111));
            }
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
	 * threads inserting the same absent k-mer may both observe it as new, so {@link #getEntries()} may
	 * marginally over-count; this never affects the filter's membership answers.
	 *
	 * @param data the k-mer, encoded as a {@code long}, to add
	 * @return {@code true} if the k-mer was not already present (at least one of its bits was newly
	 *         set), {@code false} if it was already present
	 */
	@Override
	public boolean putLongIfAbsent(final long data) {
		boolean added = false;
		for (int i = 0; i < hashes; i++) {
			// Every bit must be set even once we know the element is new, so do not short-circuit.
			if (bitVector.setAndTestWasUnset(reduce(hash(data, i)))) {
				added = true;
			}
		}
		if (added) {
			countAdded();
		}
		return added;
	}

	/**
	 * Computes the {@code i}-th hash of the given k-mer.
	 *
	 * @param data the k-mer, encoded as a {@code long}, to hash
	 * @param i    the index of the hash function to apply
	 * @return the {@code i}-th hash of the given k-mer.
	 */
	protected abstract long hash(final long data, final int i);

	/**
	 * Maps a hash value onto a valid bit index of the backing bit vector.
	 *
	 * @param v the hash value to reduce
	 * @return the bit index in {@code [0, bits)} for the given hash value.
	 */
	protected long reduce(final long v) {
		return Math.abs(v % bits);
	}

	/**
	 * Folds any k-mers counted via the concurrent insert path into {@link #entries} before the default
	 * serialization runs, so that a deserialized filter reports the correct entry count even though the
	 * concurrent counter is transient.
	 *
	 * @param out the stream the filter is written to
	 * @throws IOException if writing fails
	 */
	private void writeObject(ObjectOutputStream out) throws IOException {
		entries += drainConcurrentEntries();
		out.defaultWriteObject();
	}
}