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

import org.apache.commons.math3.distribution.NormalDistribution;
import org.metagene.genestrip.store.Database;

import it.unimi.dsi.fastutil.objects.Object2LongMap;

/**
 * Statistical estimator for the number of distinct (unique) k-mers of a taxon and for
 * how well an observed unique k-mer count matches the expectation, following the model
 * in <a href="https://arxiv.org/pdf/1602.05822.pdf">arXiv:1602.05822</a>.
 *
 * @deprecated superseded by the counting performed during matching.
 */
@Deprecated
public class UniqueKMerEstimator {
	private long totalKMers;
	private final Object2LongMap<String> storeStats;
	private final long normalizedKMersFactor;

	/**
	 * Creates an estimator over the statistics of the given database.
	 *
	 * @param storeWrapper          the database whose per-taxon k-mer statistics are used
	 * @param normalizedKMersFactor the scaling factor applied when normalizing k-mer counts
	 */
	public UniqueKMerEstimator(Database storeWrapper, long normalizedKMersFactor) {
		storeStats = storeWrapper.getStats();
		this.normalizedKMersFactor = normalizedKMersFactor;
	}

	/**
	 * Sets the total number of matched k-mers used as the normalization base.
	 *
	 * @param totalKMers the total number of matched k-mers
	 */
	public void setTotalKMers(long totalKMers) {
		this.totalKMers = totalKMers;
	}

	/**
	 * Returns the total number of matched k-mers.
	 *
	 * @return the total number of matched k-mers
	 */
	public long getTotalKMers() {
		return totalKMers;
	}
	
	/**
	 * Returns the database coverage, i.e. the ratio of observed unique k-mers to the
	 * number of k-mers stored for the taxon in the database.
	 *
	 * @param stats the per-taxon match statistics
	 * @return the database coverage ratio
	 */
	public double getDBCoverage(CountsPerTaxid stats) {
		return ((double) stats.getUniqueKMers()) / storeStats.getLong(stats.getTaxid());
	}

	/**
	 * Returns the taxon's matched k-mer count normalized by the total number of matched
	 * k-mers and the taxon's database size, scaled by the configured normalization factor.
	 *
	 * @param stats the per-taxon match statistics
	 * @return the normalized k-mer count
	 */
	public double getNormalizedKMers(CountsPerTaxid stats) {
		long totalStoredWithCounts = storeStats.getLong(null);
		return normalizedKMersFactor * ((double) stats.getKMers()) * totalStoredWithCounts / totalKMers
				/ storeStats.getLong(stats.getTaxid());
	}

	/**
	 * Returns the expected number of distinct k-mers observed given the taxon's matched
	 * k-mer count, under the model's uniform sampling assumption.
	 *
	 * @param stats the per-taxon match statistics
	 * @return the expected number of distinct observed k-mers
	 */
	public double getExpectedUniqueKMers(CountsPerTaxid stats) {
		return storeStats.getLong(stats.getTaxid()) * getAtLeastOnceUniqueKMerProb(stats);
	}

	/**
	 * Returns the probability that a given specific k-mer of the taxon is observed at
	 * least once among the matched k-mers.
	 *
	 * @param stats the per-taxon match statistics
	 * @return the probability of observing a specific k-mer at least once
	 */
	public double getAtLeastOnceUniqueKMerProb(CountsPerTaxid stats) {
		double N = storeStats.getLong(stats.getTaxid());
		double A = stats.getKMers();
		return 1 - Math.pow(1 - 1d / N, A);
	}
	
// Was not useful in practice: Removed for now:	
	/*
	 * From https://arxiv.org/pdf/1602.05822.pdf (5)
	 */
	/**
	 * Returns the variance of the number of distinct observed k-mers under the model.
	 *
	 * @param stats the per-taxon match statistics
	 * @return the variance of the number of distinct observed k-mers
	 */
	public double getUniqueKMersVariance(CountsPerTaxid stats) {
		double N = storeStats.getLong(stats.getTaxid());
		double A = stats.getKMers();

		return N * (N - 1) * Math.pow(1 - 2 / N, A) + N * Math.pow(1 - 1 / N, A) - N * N * Math.pow(1 - 1 / N, 2 * A);
	}


	/*
	 * From https://arxiv.org/pdf/1602.05822.pdf Normal distribution approximates
	 * P(k) from (1)
	 */
	/**
	 * Returns a match score between 0 and 1 for the observed unique k-mer count, computed
	 * as the normal-approximated density at the observed count divided by the density at
	 * the expected mean.
	 *
	 * @param stats the per-taxon match statistics
	 * @return a match score between 0 and 1
	 */
	public double getUniqueKMerCountMatchScore(CountsPerTaxid stats) {
		double mean = getExpectedUniqueKMers(stats);
		double sd = Math.sqrt(getUniqueKMersVariance(stats));

		if (sd == 0) {
			return 0;
		}
		NormalDistribution nd = new NormalDistribution(mean, sd);

		// We use the densitiy at X = unique kmers. We devide it by the maximum density
		// (at the mean) to obtain a score value between 0 and 1.
		return nd.density(stats.getUniqueKMers()) / nd.density(mean);

// Old attempt: P(X <= mean - |mean - unique kmers| or X >= mean + |mean - unique kmers|)
//		return nd.cumulativeProbability(mean - Math.abs(mean - stats.getUniqueKMers())) * 2;
	}

	/*
	 * From https://arxiv.org/pdf/1602.05822.pdf (9)
	 */
	/**
	 * Whether the taxon's database size and matched k-mer count lie in the range where
	 * the normal approximation of the model is considered valid.
	 *
	 * @param stats the per-taxon match statistics
	 * @return {@code true} if the normal approximation is considered valid
	 */
	public boolean isProbEstimateInRange(CountsPerTaxid stats) {
		double N = storeStats.getLong(stats.getTaxid());
		double A = stats.getKMers();

		if (N <= 5 || A <= 5) {
			return false;
		}

		return 1.4 * Math.pow(N, 0.67) <= A && A <= 1.13 * Math.pow(N, 1.19);
	}
}
