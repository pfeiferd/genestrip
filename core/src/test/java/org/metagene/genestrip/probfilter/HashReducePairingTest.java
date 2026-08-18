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

import org.junit.Test;

import java.util.HashSet;
import java.util.Random;
import java.util.Set;

import static org.junit.Assert.assertTrue;

/**
 * Pins the one invariant every filter of this package is built on: a hash and a reduction are chosen
 * together. The mixing filters carry a k-mer's low-bit entropy upwards and may therefore reduce by a
 * multiply-shift; the {@code XOR} variants do not mix and must reduce by a modulo.
 * <p>
 * The tests measure on keys with {@code 2k} significant bits, as k-mers have. That matters: with random
 * 64-bit keys a mismatched pairing looks harmless, which is how it can survive review. Here {@code k}
 * is 16, where the effect is most pronounced.
 */
public class HashReducePairingTest {
	/** Entries inserted per measurement, and lookups of absent keys used to measure the rate. */
	private static final int SIZE = 200000;
	/** Significant bases per key; 16 leaves only 32 significant bits, which exposes the pairing. */
	private static final int K = 16;
	/** Rate every correctly paired filter of this package stays below at 10 to 16 bits per key. */
	private static final double BUDGET = 3.0;

	/**
	 * Every concrete {@link ProbFilter} of this package, each on the backing its constructor picks for a
	 * filter of this size. A new implementation belongs here, or its pairing goes unmeasured.
	 */
	@Test
	public void testCorrectlyPairedVariantsStayWithinBudget() {
		assertBelow("BloomFilter", new BloomFilter(0.01, SIZE), BUDGET);
		assertBelow("MurmurBloomFilter", new MurmurBloomFilter(0.01, SIZE), BUDGET);
		assertBelow("XORBloomFilter", new XORBloomFilter(0.01, SIZE), BUDGET);
		assertBelow("BlockedBloomFilter", new BlockedBloomFilter(SIZE, 10), BUDGET);
		assertBelow("XORBlockedBloomFilter", new XORBlockedBloomFilter(SIZE, 10), BUDGET);
		assertBelow("SingleWordBloomFilter", new SingleWordBloomFilter(SIZE, 16), BUDGET);
		assertBelow("XORSingleWordBloomFilter", new XORSingleWordBloomFilter(SIZE, 16), BUDGET);
	}

	/**
	 * The bucketed backing reduces by {@code reduce} where the small one uses {@code reduceInt}, so both
	 * overrides have to be in place for the {@code XOR} variants. Only the two one-word families appear
	 * here: the {@link BloomFilter} family has a single backing, the
	 * {@link org.metagene.genestrip.util.LargeBitVector} it always uses, so the sizings above already
	 * exercise the only reduction it has.
	 */
	@Test
	public void testCorrectlyPairedVariantsStayWithinBudgetOnLargeBacking() {
		assertBelow("BlockedBloomFilter", BlockedBloomFilter.newLargeBacked(SIZE, 10), BUDGET);
		assertBelow("XORBlockedBloomFilter", XORBlockedBloomFilter.newLargeBackedXOR(SIZE, 10), BUDGET);
		assertBelow("SingleWordBloomFilter", SingleWordBloomFilter.newLargeBacked(SIZE, 16), BUDGET);
		assertBelow("XORSingleWordBloomFilter", XORSingleWordBloomFilter.newLargeBackedXOR(SIZE, 16), BUDGET);
	}

	/**
	 * An XOR hash left on the inherited multiply-shift must measure clearly worse than the same filter
	 * with the modulo its variant uses, in every family. Were this to stop holding, the reason for the
	 * {@code XOR} subclasses would have gone with it.
	 */
	@Test
	public void testMismatchedPairingIsWorse() {
		assertWorse("BloomFilter", new Mismatched(0.01, SIZE), new XORBloomFilter(0.01, SIZE));
		assertWorse("BlockedBloomFilter", new MismatchedBlocked(SIZE, 10), new XORBlockedBloomFilter(SIZE, 10));
		assertWorse("SingleWordBloomFilter", new MismatchedSingleWord(SIZE, 16),
				new XORSingleWordBloomFilter(SIZE, 16));
	}

	private void assertBelow(String name, ProbFilter filter, double budget) {
		double rate = measure(filter);
		assertTrue(name + " false positive rate: " + rate + "%", rate < budget);
	}

	private void assertWorse(String family, ProbFilter mismatched, ProbFilter paired) {
		double bad = measure(mismatched);
		double good = measure(paired);
		assertTrue(family + ": xor hash on the multiply-shift measured " + bad + "%, with the modulo " + good
				+ "% - the modulo is supposed to be clearly better", bad > 1.5 * good);
	}

	/**
	 * Inserts {@link #SIZE} k-mer-like keys and returns the percentage of {@link #SIZE} lookups of keys
	 * that were never inserted which the filter nevertheless reports as present.
	 *
	 * @param filter the filter to measure
	 * @return the measured false-positive rate in percent
	 */
	private double measure(ProbFilter filter) {
		Random random = new Random(2026);
		long mask = (1L << (2 * K)) - 1;
		Set<Long> inserted = new HashSet<>();
		while (inserted.size() < SIZE) {
			inserted.add(random.nextLong() & mask);
		}
		for (long key : inserted) {
			filter.putLong(key);
		}
		int hits = 0;
		int lookups = 0;
		while (lookups < SIZE) {
			long key = random.nextLong() & mask;
			if (inserted.contains(key)) {
				continue;
			}
			lookups++;
			if (filter.containsLong(key)) {
				hits++;
			}
		}
		return 100.0 * hits / lookups;
	}

	/** {@link XORBloomFilter}'s hash without its modulo, i.e. the pairing the design forbids. */
	private static class Mismatched extends BloomFilter {
		private static final long serialVersionUID = 1L;

		Mismatched(double fpp, long expectedInsertions) {
			super(fpp, expectedInsertions);
		}

		@Override
		protected long hash(long data, int i) {
			return hashFactors[i] ^ data;
		}
	}

	/** {@link XORBlockedBloomFilter}'s hash without its modulo. */
	private static class MismatchedBlocked extends BlockedBloomFilter {
		private static final long serialVersionUID = 1L;

		MismatchedBlocked(long expectedInsertions, int bitsPerKey) {
			super(expectedInsertions, bitsPerKey);
		}

		@Override
		protected long hash(long x) {
			return seed ^ x;
		}
	}

	/** {@link XORSingleWordBloomFilter}'s hash without its modulo. */
	private static class MismatchedSingleWord extends SingleWordBloomFilter {
		private static final long serialVersionUID = 1L;

		MismatchedSingleWord(long expectedInsertions, int bitsPerKey) {
			super(expectedInsertions, bitsPerKey);
		}

		@Override
		protected long hash(long x) {
			return seed ^ x;
		}
	}
}
