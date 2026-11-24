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

import java.util.Random;

import org.junit.Test;
import org.metagene.genestrip.match.SharedFractionSimulator.Counts;
import org.metagene.genestrip.match.SharedFractionSimulator.Isolate;

/**
 * Checks {@link SharedFractionEstimator} against simulated isolates (see {@link SharedFractionSimulator}).
 * The claims tested are those of {@code SharedFraction.md}: both estimators recover the shared
 * fraction under uniform coverage, their intervals hold their level only when the runs of hits are
 * accounted for, random sequencing errors do not bias them, and uneven coverage biases both downwards,
 * Chao's estimator less - and a multi-copy plasmid only the one-bit estimator.
 */
public class SharedFractionSimulationTest {
	private static final int G = 300_000;

	private final SharedFractionEstimator estimator = new SharedFractionEstimator(SharedFractionSimulator.K, 0.95);

	// The node and the isolate are drawn from one random sequence, the isolate with 8 % of the genome
	// deleted and one SNP per 1000 bases.
	private static Isolate isolate(long seed) {
		Random random = new Random(seed);
		return new SharedFractionSimulator(random, G, 0.3, 2000).isolate(random, 0.08, 5000, 1e-3);
	}

	@Test
	public void testUniformCoverageRecoversPi() {
		Isolate isolate = isolate(1);
		for (double x : new double[] { 1, 2, 4 }) {
			Counts counts = isolate.sample(new Random(2), x, null, 0);
			assertWithin("one-bit, x = " + x, isolate.pi, counts.endpoint(estimator));
			assertWithin("Chao, x = " + x, isolate.pi, counts.chao(estimator));
		}
	}

	@Test
	public void testRandomErrorsDoNotBias() {
		Isolate isolate = isolate(3);
		Counts counts = isolate.sample(new Random(4), 2, null, 0.01);
		assertWithin("one-bit", isolate.pi, counts.endpoint(estimator));
		assertWithin("Chao", isolate.pi, counts.chao(estimator));
	}

	// The 95 % intervals with b over runs must contain pi in about 95 % of the samples; with the hits
	// taken for independent draws (b = 1) they are far too narrow.
	@Test
	public void testIntervalCoverage() {
		Isolate isolate = isolate(5);
		Random random = new Random(6);
		int samples = 200;
		int endpoint = 0, chao = 0, chaoIndependent = 0;
		for (int i = 0; i < samples; i++) {
			Counts counts = isolate.sample(random, 1, null, 0);
			endpoint += covers(counts.endpoint(estimator), isolate.pi);
			chao += covers(counts.chao(estimator), isolate.pi);
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
		Isolate isolate = isolate(7);
		Counts counts = isolate.sample(new Random(8), 1, isolate.plasmidWeights(new Random(9), 0.05, 20), 0);
		double endpointError = counts.endpoint(estimator)[0] - isolate.pi;
		assertTrue("one-bit error " + endpointError, endpointError < -0.2);
		assertWithin("Chao", isolate.pi, counts.chao(estimator));
	}

	// Read-start density constant on 1 kb windows and gamma-distributed across them (shape 2, CV 0.71).
	@Test
	public void testUnevenCoverage() {
		Isolate isolate = isolate(10);
		Counts counts = isolate.sample(new Random(11), 2, isolate.gammaWeights(new Random(12), 2), 0);
		double endpointError = counts.endpoint(estimator)[0] - isolate.pi;
		double chaoError = counts.chao(estimator)[0] - isolate.pi;
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
}
