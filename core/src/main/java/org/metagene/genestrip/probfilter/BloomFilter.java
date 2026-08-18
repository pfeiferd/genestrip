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
package org.metagene.genestrip.probfilter;

import org.metagene.genestrip.util.LargeBitVector;

import java.util.Random;

/**
 * The classical {@link ProbFilter}, backed by a {@link LargeBitVector}: it sizes the bit vector and
 * the number of hash functions from the expected insertions and target false-positive probability, and
 * sets one bit per hash function. Unlike {@link BlockedBloomFilter} and
 * {@link SingleWordBloomFilter}, which confine a key to one or two adjacent words, those bits are
 * spread over the whole vector, which buys the best rate per bit and costs one memory access per hash
 * function.
 * <p>
 * This class is the <em>mixing</em> variant of its family, i.e. the counterpart of
 * {@link XORBloomFilter}. Like every filter of this package it exists in two variants that differ
 * in one pair of choices which must always be made together:
 * <ul>
 * <li>a mixing hash - {@link #hash(long, int)} here, MurmurHash3 in {@link MurmurBloomFilter} -
 * which spreads a k-mer's entropy over the whole word and may therefore use the multiply-shift
 * {@link #reduce(long)};</li>
 * <li>a bare exclusive or - {@link XORBloomFilter} - which must therefore reduce by a modulo, as it
 * overrides {@link #reduce(long)} to do.</li>
 * </ul>
 * Pairing a non-mixing hash with a multiply-shift reduction costs orders of magnitude of
 * false-positive rate, so the two always travel together; see {@link XORBloomFilter#reduce(long)}
 * for the measurements.
 */
public class BloomFilter implements ProbFilter {
	private static final long serialVersionUID = 2L;

	/** Number of Simpson steps used by {@link #estimateDistinctValues(long)}; even, as the rule requires. */
	private static final int ESTIMATION_STEPS = 2000;
	/** Cap on the rate used there, so that a saturated filter yields a large value rather than infinity. */
	private static final double MAX_ESTIMATION_FPP = 0.999999;

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
	public BloomFilter(double fpp, long expectedInsertions) {
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
	 * Returns the false-positive probability the filter has reached after the given number of
	 * insertions.
	 * <p>
	 * This is not the same as {@link #getFpp()}, which is the probability the filter was
	 * <em>sized</em> for: that one is reached at {@link #getExpectedInsertions()} insertions, stays
	 * below it for fewer and rises above it for more. A filter that is filled beyond what it was
	 * built for still answers, it merely answers wrongly more often, and this method is what says by
	 * how much.
	 * <p>
	 * The count has to be supplied because the filter cannot know it. {@link #putLong(long)} does
	 * report whether a value was new, but nothing keeps a tally, and the caller is in any case the
	 * only one that can say how many <em>distinct</em> values were offered - inserting the same
	 * value twice sets the same bits and does not raise the probability at all.
	 *
	 * @param insertions the number of distinct values inserted so far
	 * @return the probability that {@link #containsLong(long)} answers {@code true} for a value that
	 *         was never inserted
	 * @throws IllegalArgumentException if {@code insertions} is negative
	 */
	public double getFpp(long insertions) {
		if (insertions < 0) {
			throw new IllegalArgumentException("insertions must be >= 0");
		}
		// The expected share of bits that the k hashes of n values have set, and then the probability
		// that the k hashes of an absent value all land on set bits. expm1 rather than 1 - exp: for a
		// filter far from full the exponent is close to zero, and there the difference of the two is
		// all that is left of the result.
		double setRatio = -Math.expm1(-((double) hashes * insertions) / bits);
		return Math.pow(setRatio, hashes);
	}

	/**
	 * Returns how many distinct values are estimated to have been offered to this filter, given that
	 * it reported {@code counted} of them as new.
	 * <p>
	 * This is the counterpart of {@link #getFpp(long)}: that one goes from a number of insertions to
	 * a false-positive rate, this one from an observed count back to the number of values behind it.
	 * A caller that counts how often {@link #putLong(long)} reported a value as new undercounts,
	 * because a value the filter falsely holds to be present is never counted, and how badly depends
	 * on how full the filter was at the time - not on what it was sized for, and not on where it
	 * ended up.
	 * <p>
	 * With {@code c} counted so far, a further value is counted unless the filter hides it, so
	 * {@code dc/dT = 1 - fpp(c)} and the number truly offered is the integral of
	 * {@code 1 / (1 - fpp(c))} over {@code c}, evaluated here by Simpson's rule. For a filter well
	 * inside its sizing this differs little from dividing by {@code 1 - fpp}; beyond it, the two part
	 * company sharply - at four times the sizing the constant-rate correction is out by some 14 %
	 * where this one is out by 4 %.
	 * <p>
	 * <strong>How far it can be trusted.</strong> The estimate is not equally good everywhere, and
	 * where it fails it does so in one direction: it <em>under</em>-estimates, because a saturated
	 * filter stops leaving a trace of new values at all and no arithmetic recovers what was never
	 * recorded. Measured against filters filled with a known number of distinct values, sized for
	 * 200,000 at a target rate of 1 %:
	 * <ul>
	 * <li>up to the sizing: within 0.01 %;</li>
	 * <li>at twice the sizing: within 0.3 %;</li>
	 * <li>at three times: about 3 % low;</li>
	 * <li>at four times: about 9 % low;</li>
	 * <li>at eight times: about 43 % low, i.e. no longer an estimate.</li>
	 * </ul>
	 * A filter built for a lower rate holds up longer - at a target of 0.1 % the same points are
	 * 0.00 %, 0.00 %, 0.7 % and 4 % - because it has more bits and more hash functions to saturate.
	 * What decides the accuracy is not the multiple of the sizing as such but how full the filter
	 * ended up: while {@link #getFpp(long)} of the counted values stays below roughly 0.15 the
	 * estimate is good to a fraction of a percent, and it degrades from there. A caller that cares
	 * should test that rate and treat the result as a lower bound above it.
	 *
	 * @param counted the number of values this filter reported as new
	 * @return the estimated number of distinct values offered, never less than {@code counted}
	 * @throws IllegalArgumentException if {@code counted} is negative
	 */
	public double estimateDistinctValues(long counted) {
		if (counted < 0) {
			throw new IllegalArgumentException("counted must be >= 0");
		}
		if (counted == 0) {
			return 0;
		}
		double h = (double) counted / ESTIMATION_STEPS;
		double sum = 0;
		for (int i = 0; i <= ESTIMATION_STEPS; i++) {
			// A saturated filter would divide by zero here; the cap turns that into a large but
			// finite value. It is reached only once the filter hides values faster than any
			// correction recovers them, where the result is a lower bound whatever it is.
			double value = 1d / (1d - Math.min(getFpp((long) (i * h)), MAX_ESTIMATION_FPP));
			sum += (i == 0 || i == ESTIMATION_STEPS) ? value : (i % 2 == 1 ? 4 * value : 2 * value);
		}
		return Math.max(counted, sum * h / 3);
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
	 * Computes the {@code i}-th hash of the given k-mer by a bit-mixing hash (moved here from the former
	 * {@code LemireOptBloomFilter}) that carries a k-mer's low-bit entropy upwards and hence pairs with
	 * the multiply-shift {@link #reduce(long)}. A subclass overriding it with a hash that does not mix
	 * must override {@link #reduce(long)} with a modulo as well - see {@link XORBloomFilter}.
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
	 * Maps a hash value onto a valid bit index of the backing bit vector, using Lemire's multiply-shift
	 * alternative to the modulo widened to 64 bits: the upper half of the 128-bit product of the hash
	 * (read as unsigned) with {@link #bits}, which {@link Math#multiplyHigh(long, long)} yields, lands
	 * uniformly in {@code [0, bits)}. That avoids the 64-bit division a modulo costs per hash function.
	 * <p>
	 * This is the reduction of the <em>mixing</em> filters of this package, and it is only sound for a
	 * hash that mixes: every multiply-shift is driven by the <em>high</em> bits of its input, whereas a
	 * k-mer carries its entropy in the low ones. {@link #hash(long, int)} and the hash of
	 * {@link MurmurBloomFilter} both carry that entropy upwards, so they may reduce this way.
	 * A hash that does not mix must not - see {@link XORBloomFilter#reduce(long)}, which overrides
	 * this with a modulo for exactly that reason.
	 *
	 * @param v the hash value to reduce
	 * @return the bit index in {@code [0, bits)} for the given hash value.
	 */
	protected long reduce(final long v) {
		// Math.multiplyHigh is signed, so a negative left operand needs the range added back to reach
		// the unsigned product's upper half; 'bits' is always positive, so only that side corrects.
		return Math.multiplyHigh(v, bits) + ((v >> 63) & bits);
	}
}
