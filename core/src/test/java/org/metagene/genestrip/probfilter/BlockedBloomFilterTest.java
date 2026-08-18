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

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import org.junit.Test;

public class BlockedBloomFilterTest extends BloomFilterTest {
	@Override
	protected double createFpp() {
		return 0.01;
	}

	@Override
	protected ProbFilter createFilter(long size, double fpp) {
		// The filter picks its backing (small array vs. bucketed grid) from the size relative to
		// MAX_SMALL_CAPACITY at construction, which no test-sized filter reaches, so ask for the
		// bucketed backing explicitly.
		return isTestLarge() ? BlockedBloomFilter.newLargeBacked(size, 10)
				: new BlockedBloomFilter(size, 10);
	}

	@Override
	protected boolean isTestLarge() {
		return false;
	}

	/** Mixes consecutive numbers over the whole word, so the filter sees no structure. */
	private static long mixed(long i) {
		long z = i * 0x9E3779B97F4A7C15L;
		z = (z ^ (z >>> 30)) * 0xBF58476D1CE4E5B9L;
		return z ^ (z >>> 27);
	}

	/**
	 * The rate the filter reports for a number of insertions must be the rate it actually shows once
	 * that many are in it. The tolerance is generous below the sizing, where the formula is known to
	 * underestimate and both figures are negligible, and tight from the sizing onwards, which is where
	 * a caller decides anything on it.
	 */
	@Test
	public void testReportedFppMatchesAFilledFilter() {
		long sizedFor = 200_000L;
		for (double load : new double[] { 1.0, 2.0, 4.0 }) {
			long insertions = (long) (load * sizedFor);
			BlockedBloomFilter filter = new BlockedBloomFilter(sizedFor);
			for (long i = 0; i < insertions; i++) {
				filter.putLong(mixed(i));
			}
			int probes = 200_000, hits = 0;
			for (long i = 1; i <= probes; i++) {
				if (filter.containsLong(mixed(-i))) {
					hits++;
				}
			}
			double measured = (double) hits / probes;
			double reported = filter.getFpp(insertions);
			assertEquals("at " + load + " times the sizing", measured, reported, 0.005 + 0.2 * measured);
		}
	}

	/** An empty filter answers no membership, more insertions never lower the rate, and it stays a probability. */
	@Test
	public void testFppGrowsWithInsertionsAndStaysAProbability() {
		BlockedBloomFilter filter = new BlockedBloomFilter(100_000L);
		assertEquals(0.0, filter.getFpp(0), 0.0);
		double previous = -1;
		for (long insertions = 0; insertions <= 1_000_000L; insertions += 50_000L) {
			double current = filter.getFpp(insertions);
			assertTrue("the rate fell at " + insertions, current >= previous);
			assertTrue("the rate left [0,1] at " + insertions, current >= 0 && current <= 1);
			previous = current;
		}
		try {
			filter.getFpp(-1);
			fail("a negative number of insertions should have been rejected");
		} catch (IllegalArgumentException expected) {
			// as it should be
		}
	}

	/**
	 * Pins the margin the sizing check of {@code kmerindexbloom} rests on. That goal builds the filter
	 * for an estimated number of entries and refuses it once its false-positive rate has grown to 1.75
	 * times the rate it would have at that number - which has to be far above what the estimate's own
	 * error can produce and far below what a filter sized by a genuinely wrong number reaches.
	 */
	@Test
	public void testTheSizingMarginSeparatesNoiseFromRealMisses() {
		final double factor = 1.75;
		final long sized = 100_000_000L;
		BlockedBloomFilter filter = new BlockedBloomFilter(sized);
		double designFpp = filter.getFpp(sized);
		double gate = designFpp * factor;

		// The estimate is good to well under a per cent, so these must pass.
		for (double tooLow : new double[] { 0.002, 0.005, 0.01, 0.05 }) {
			assertTrue("an estimate " + (100 * tooLow) + "% short should not be refused (fpp "
					+ filter.getFpp((long) (sized * (1 + tooLow))) + " against a gate of " + gate + ")",
					filter.getFpp((long) (sized * (1 + tooLow))) <= gate);
		}
		// An estimate that is wrong rather than noisy is wrong by factors, and must be caught.
		for (double tooLow : new double[] { 1.0, 3.0, 10.0 }) {
			assertTrue("an estimate " + (100 * tooLow) + "% short should be refused",
					filter.getFpp((long) (sized * (1 + tooLow))) > gate);
		}
		// And the gate must not fire on a filter that got exactly what it was built for.
		assertTrue("a filter holding exactly its sizing count must not be refused",
				filter.getFpp(sized) <= gate);
	}
}
