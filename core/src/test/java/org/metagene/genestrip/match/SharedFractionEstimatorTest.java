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

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

public class SharedFractionEstimatorTest {
	@Test
	public void testNormalQuantile() {
		assertEquals(0, SharedFractionEstimator.normalQuantile(0.5), 1e-9);
		assertEquals(1.959964, SharedFractionEstimator.normalQuantile(0.975), 1e-6);
		assertEquals(-1.959964, SharedFractionEstimator.normalQuantile(0.025), 1e-6);
		assertEquals(2.575829, SharedFractionEstimator.normalQuantile(0.995), 1e-6);
		assertEquals(-3.090232, SharedFractionEstimator.normalQuantile(0.001), 1e-6);
	}

	// With the counts at their expectations under equal rates, both estimators recover M / N.
	@Test
	public void testBothEstimatorsAtExpectation() {
		long n = 1_000_000;
		long m = 900_000;
		for (double lambda : new double[] { 0.5, 1, 2, 4 }) {
			long c = Math.round(lambda * m);
			long u = Math.round(m * (1 - Math.exp(-lambda)));
			long f1 = Math.round(m * lambda * Math.exp(-lambda));
			long f2 = Math.round(m * lambda * lambda * Math.exp(-lambda) / 2);
			assertEquals("endpoint, lambda " + lambda, 0.9, SharedFractionEstimator.endpoint(c, u, n), 1e-4);
			assertEquals("chao, lambda " + lambda, 0.9, SharedFractionEstimator.chao(u, f1, f2) / n, 1e-3);
		}
	}

	// For b = 1 the variance of the Chao estimate reduces to Chao's (1987) variance for independent draws,
	// f2 (r^4 / 4 + r^3 + r^2 / 2) with r = f1 / f2.
	@Test
	public void testChaoVarianceForIndependentDraws() {
		long n = 1_000_000;
		long m = 900_000;
		for (double lambda : new double[] { 0.5, 1, 2 }) {
			long u = Math.round(m * (1 - Math.exp(-lambda)));
			long f1 = Math.round(m * lambda * Math.exp(-lambda));
			long f2 = Math.round(m * lambda * lambda * Math.exp(-lambda) / 2);
			double r = (double) f1 / f2;
			double chao1987 = f2 * (r * r * r * r / 4 + r * r * r + r * r / 2);
			double se = SharedFractionEstimator.chaoStandardError(u, f1, f2, n, 1);
			assertEquals("lambda " + lambda, 1, se * se * n * n / chao1987, 0.01);
		}
	}

	// Runs of b hits widen both standard errors by roughly sqrt(b).
	@Test
	public void testRunsWidenTheStandardErrors() {
		long n = 1_000_000;
		long m = 900_000;
		double lambda = 1;
		long c = Math.round(lambda * m);
		long u = Math.round(m * (1 - Math.exp(-lambda)));
		long f1 = Math.round(m * lambda * Math.exp(-lambda));
		long f2 = Math.round(m * lambda * lambda * Math.exp(-lambda) / 2);
		double ratioEndpoint = SharedFractionEstimator.endpointStandardError(c, u, n, 100)
				/ SharedFractionEstimator.endpointStandardError(c, u, n, 1);
		double ratioChao = SharedFractionEstimator.chaoStandardError(u, f1, f2, n, 100)
				/ SharedFractionEstimator.chaoStandardError(u, f1, f2, n, 1);
		assertTrue("endpoint ratio " + ratioEndpoint, ratioEndpoint > 4 && ratioEndpoint < 12);
		assertTrue("chao ratio " + ratioChao, ratioChao > 4 && ratioChao < 12);
	}

	// No overflow at high depth: the standard errors stay finite and tiny.
	@Test
	public void testHighDepth() {
		long n = 1_000_000;
		long m = 900_000;
		long c = 500L * m;
		assertEquals(0.9, SharedFractionEstimator.endpoint(c, m, n), 1e-9);
		double se = SharedFractionEstimator.endpointStandardError(c, m - 1, n, 120);
		assertTrue("se " + se, se >= 0 && se < 1e-6);
		double chaoSe = SharedFractionEstimator.chaoStandardError(m - 2, 1, 1, n, 120);
		assertTrue("chao se " + chaoSe, chaoSe >= 0 && chaoSe < 1e-3);
	}

	// estimate() fills the counts: the interval brackets the estimate, stays within [u/N, 1], and the
	// distance interval is the image of the fraction interval.
	@Test
	public void testEstimateFillsCounts() {
		for (boolean counted : new boolean[] { false, true }) {
			CountsPerTaxid counts = new CountsPerTaxid(0, "1", 0);
			counts.kmers = 900_000;
			counts.uniqueKmers = 560_000;
			counts.contigs = 9_000;
			counts.contigLenSquaredSum = 9_000L * 100 * 100;
			if (counted) {
				counts.kmersHitOnce = 330_000;
				counts.kmersHitTwice = 165_000;
			}
			counts.completeValues(1, 1_000_000, null);
			SharedFractionEstimator estimator = new SharedFractionEstimator(31, 0.95);
			estimator.estimate(counts);
			double pi = counts.getSharedFraction();
			assertTrue(counts.getSharedFractionLow() <= pi && pi <= counts.getSharedFractionHigh());
			assertTrue(counts.getSharedFractionLow() >= 0.56 && counts.getSharedFractionHigh() <= 1);
			assertEquals(1 - Math.pow(pi, 1.0 / 31), counts.getEvoDistance(), 1e-12);
			assertEquals(1 - Math.pow(counts.getSharedFractionHigh(), 1.0 / 31), counts.getEvoDistanceLow(), 1e-12);
			assertEquals(1 - Math.pow(counts.getSharedFractionLow(), 1.0 / 31), counts.getEvoDistanceHigh(), 1e-12);
			assertEquals(Math.round((1 - pi) * 1_000_000), (long) counts.getMissingKMers());
		}
	}
}
