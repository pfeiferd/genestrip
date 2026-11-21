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
package org.metagene.genestrip.match;

import static org.junit.Assert.assertTrue;

import java.util.Arrays;
import java.util.Random;

import org.junit.Test;

/**
 * Checks {@link SharedFractionEstimator} against simulated isolates, as described in
 * {@code SharedFraction.md}: a circular genome carries runs of k-mers specific to a database node, the
 * isolate lacks some of them through deletions and SNPs, and reads drawn from it hit the node's k-mers
 * in runs. The claims tested are those of the documentation: both estimators recover the shared
 * fraction under uniform coverage, their intervals hold their level only when the runs of hits are
 * accounted for, random sequencing errors do not bias them, and uneven coverage biases both downwards,
 * Chao's estimator less - and a multi-copy plasmid only the one-bit estimator.
 */
public class SharedFractionSimulationTest {
	private static final int G = 300_000;
	private static final int K = 31;
	private static final int L = 150;
	private static final int W = L - K + 1;

	private final SharedFractionEstimator estimator = new SharedFractionEstimator(K, 0.95);

	@Test
	public void testUniformCoverageRecoversPi() {
		Isolate isolate = new Isolate(1);
		for (double x : new double[] { 1, 2, 4 }) {
			Counts counts = isolate.sample(new Random(2), x, null, 0);
			assertWithin("one-bit, x = " + x, isolate.pi, counts.endpoint());
			assertWithin("Chao, x = " + x, isolate.pi, counts.chao());
		}
	}

	@Test
	public void testRandomErrorsDoNotBias() {
		Isolate isolate = new Isolate(3);
		Counts counts = isolate.sample(new Random(4), 2, null, 0.01);
		assertWithin("one-bit", isolate.pi, counts.endpoint());
		assertWithin("Chao", isolate.pi, counts.chao());
	}

	// The 95 % intervals with b over runs must contain pi in about 95 % of the samples; with the hits
	// taken for independent draws (b = 1) they are far too narrow.
	@Test
	public void testIntervalCoverage() {
		Isolate isolate = new Isolate(5);
		Random random = new Random(6);
		int samples = 200;
		int endpoint = 0, chao = 0, chaoIndependent = 0;
		for (int i = 0; i < samples; i++) {
			Counts counts = isolate.sample(random, 1, null, 0);
			endpoint += covers(counts.endpoint(), isolate.pi);
			chao += covers(counts.chao(), isolate.pi);
			chaoIndependent += covers(estimator.chaoInterval(counts.u, counts.f1, counts.f2, counts.n, 1),
					isolate.pi);
		}
		assertTrue("one-bit coverage " + endpoint, endpoint >= 0.85 * samples && endpoint <= 0.995 * samples);
		assertTrue("Chao coverage " + chao, chao >= 0.85 * samples && chao <= 0.995 * samples);
		assertTrue("Chao coverage with b = 1: " + chaoIndependent, chaoIndependent < 0.5 * samples);
	}

	// 5 % of the shared k-mers on a plasmid of 20 copies.
	@Test
	public void testPlasmid() {
		Isolate isolate = new Isolate(7);
		Counts counts = isolate.sample(new Random(8), 1, isolate.plasmidWeights(new Random(9), 0.05, 20), 0);
		assertTrue("one-bit error " + (counts.endpoint()[0] - isolate.pi), counts.endpoint()[0] - isolate.pi < -0.2);
		assertWithin("Chao", isolate.pi, counts.chao());
	}

	// Read-start density constant on 1 kb windows and gamma-distributed across them (shape 2, CV 0.71).
	@Test
	public void testUnevenCoverage() {
		Isolate isolate = new Isolate(10);
		Counts counts = isolate.sample(new Random(11), 2, isolate.gammaWeights(new Random(12), 2), 0);
		double endpointError = counts.endpoint()[0] - isolate.pi;
		double chaoError = counts.chao()[0] - isolate.pi;
		assertTrue("one-bit error " + endpointError, endpointError < -0.03);
		assertTrue("Chao error " + chaoError, chaoError < 0 && chaoError > endpointError / 1.5);
	}

	// The estimate must lie within three of its own standard errors of pi; the standard error is read off
	// the 95 % interval, whose half-width is 1.96 of them (where the interval is not cut at u/N or 1).
	private static void assertWithin(String message, double pi, double[] estimate) {
		double se = (estimate[2] - estimate[1]) / (2 * 1.959964);
		assertTrue(message + ": estimate " + estimate[0] + ", pi " + pi + ", se " + se,
				Math.abs(estimate[0] - pi) <= 3 * se);
	}

	private static int covers(double[] interval, double pi) {
		return interval[1] <= pi && pi <= interval[2] ? 1 : 0;
	}

	/** The counts after a sample of reads, as the matcher would report them for the node. */
	private class Counts {
		long n, c, u, f1, f2, runs, runSquares;

		double[] endpoint() {
			return estimator.endpointInterval(c, u, n, SharedFractionEstimator.runLength(runSquares, runs));
		}

		double[] chao() {
			return estimator.chaoInterval(u, f1, f2, n, SharedFractionEstimator.runLength(runSquares, runs));
		}
	}

	/** A genome with the node's k-mers, and the isolate's copy of it. */
	private class Isolate {
		final boolean[] specific = new boolean[G];
		final boolean[] deleted = new boolean[G];
		final boolean[] present = new boolean[G];
		final int n;
		final double pi;

		Isolate(long seed) {
			Random random = new Random(seed);
			runs(random, specific, 0.3, 2000);
			runs(random, deleted, 0.08, 5000);
			boolean[] destroyed = new boolean[G];
			for (int i = 0; i < G; i++) {
				if (random.nextDouble() < 1e-3) {
					// A SNP at base i destroys the k-mers starting at i - k + 1 .. i.
					for (int d = 0; d < K; d++) {
						destroyed[Math.floorMod(i - d, G)] = true;
					}
				}
			}
			int count = 0, shared = 0;
			for (int i = 0; i < G; i++) {
				present[i] = specific[i] && !deleted[i] && !destroyed[i];
				count += specific[i] ? 1 : 0;
				shared += present[i] ? 1 : 0;
			}
			n = count;
			pi = (double) shared / count;
		}

		private void runs(Random random, boolean[] mask, double share, int meanLength) {
			int marked = 0;
			while (marked < share * G) {
				int start = random.nextInt(G);
				int length = Math.max(1, (int) (-meanLength * Math.log(1 - random.nextDouble())));
				for (int i = 0; i < length; i++) {
					int p = (start + i) % G;
					if (!mask[p]) {
						mask[p] = true;
						marked++;
					}
				}
			}
		}

		double[] plasmidWeights(Random random, double share, double copies) {
			double[] weights = new double[G];
			Arrays.fill(weights, 1);
			int[] shared = new int[G];
			int m = 0;
			for (int i = 0; i < G; i++) {
				if (present[i]) {
					shared[m++] = i;
				}
			}
			int first = random.nextInt(m - (int) (share * m));
			for (int i = shared[first]; i <= shared[first + (int) (share * m)]; i++) {
				weights[i] = copies;
			}
			return weights;
		}

		double[] gammaWeights(Random random, double shape) {
			double[] weights = new double[G];
			for (int window = 0; window * 1000 < G; window++) {
				double w = gamma(random, shape) / shape;
				for (int i = window * 1000; i < Math.min(G, (window + 1) * 1000); i++) {
					weights[i] = w;
				}
			}
			return weights;
		}

		/**
		 * Draws reads until the hits reach {@code x} per stored k-mer. Reads start at random,
		 * proportionally to {@code weights} if given; reads from a deleted segment do not exist. With an
		 * error rate, each base is wrong with that probability, and a k-mer spanning a wrong base does not
		 * hit.
		 */
		Counts sample(Random random, double x, double[] weights, double errorRate) {
			double[] cumulative = new double[G];
			double total = 0;
			for (int i = 0; i < G; i++) {
				total += deleted[i] ? 0 : (weights == null ? 1 : weights[i]);
				cumulative[i] = total;
			}
			int[] hits = new int[G];
			boolean[] error = new boolean[L];
			Counts counts = new Counts();
			counts.n = n;
			while (counts.c < x * n) {
				int start = Arrays.binarySearch(cumulative, random.nextDouble() * total);
				start = start < 0 ? -start - 1 : start;
				if (deleted[start] || deleted[(start + L - 1) % G]) {
					continue;
				}
				for (int j = 0; j < L; j++) {
					error[j] = errorRate > 0 && random.nextDouble() < errorRate;
				}
				int errorsInKMer = 0;
				for (int j = 0; j < K - 1; j++) {
					errorsInKMer += error[j] ? 1 : 0;
				}
				int run = 0;
				for (int i = 0; i < W; i++) {
					errorsInKMer += error[i + K - 1] ? 1 : 0;
					int p = (start + i) % G;
					if (present[p] && errorsInKMer == 0) {
						int old = hits[p]++;
						counts.c++;
						if (old == 0) {
							counts.u++;
							counts.f1++;
						} else if (old == 1) {
							counts.f1--;
							counts.f2++;
						} else if (old == 2) {
							counts.f2--;
						}
						run++;
					} else if (run > 0) {
						counts.runs += run;
						counts.runSquares += (long) run * run;
						run = 0;
					}
					errorsInKMer -= error[i] ? 1 : 0;
				}
				if (run > 0) {
					counts.runs += run;
					counts.runSquares += (long) run * run;
				}
			}
			return counts;
		}
	}

	/** Marsaglia and Tsang's method for a gamma variable with the given shape and scale one. */
	private static double gamma(Random random, double shape) {
		double d = shape - 1.0 / 3;
		double c = 1 / Math.sqrt(9 * d);
		while (true) {
			double z = random.nextGaussian();
			double v = Math.pow(1 + c * z, 3);
			if (v > 0 && Math.log(random.nextDouble()) < z * z / 2 + d - d * v + d * Math.log(v)) {
				return d * v;
			}
		}
	}
}
