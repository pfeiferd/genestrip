package org.metagene.genestrip.match;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.distribution.NormalDistribution;
import org.metagene.genestrip.match.FastqKMerMatcher.StatsPerTaxid;
import org.metagene.genestrip.store.KMerStoreWrapper.StoreStatsPerTaxid;

public class UniqueKMerEstimator {
	private long totalKMers;
	private final Map<String, StoreStatsPerTaxid> taxidToStoreStats;

	public UniqueKMerEstimator(List<StoreStatsPerTaxid> storeStatsPerTaxid) {
		taxidToStoreStats = new HashMap<String, StoreStatsPerTaxid>();
		for (StoreStatsPerTaxid stats : storeStatsPerTaxid) {
			taxidToStoreStats.put(stats.getTaxid(), stats);
		}
	}

	public void setTotalKMers(long totalKMers) {
		this.totalKMers = totalKMers;
	}

	public long getTotalKMers() {
		return totalKMers;
	}

	public double getNormalizedKMers(StatsPerTaxid stats) {
		return ((double) stats.getKMers()) / totalKMers
				/ taxidToStoreStats.get(stats.getTaxid()).assignedKMers;
	}

	/*
	 * From https://arxiv.org/pdf/1602.05822.pdf (5)
	 */
	public double getExpectedUniqueKMers(StatsPerTaxid stats) {
		return taxidToStoreStats.get(stats.getTaxid()).assignedKMers * getAtLeastOnceUniqueKMerProb(stats);
	}

	/*
	 * From https://arxiv.org/pdf/1602.05822.pdf (5)
	 */
	public double getUniqueKMersVariance(StatsPerTaxid stats) {
		double N = taxidToStoreStats.get(stats.getTaxid()).assignedKMers;
		double A = stats.getKMers();

		return N * (N - 1) * Math.pow(1 - 2 / N, A) + N * Math.pow(1 - 1 / N, A) - N * N * Math.pow(1 - 1 / N, 2 * A);
	}

	public double getAtLeastOnceUniqueKMerProb(StatsPerTaxid stats) {
		double N = taxidToStoreStats.get(stats.getTaxid()).assignedKMers;
		double A = stats.getKMers();
		return 1 - Math.pow(1 - 1d / N, A);
	}

	/*
	 * From https://arxiv.org/pdf/1602.05822.pdf Normal distribution approximates P(k) from (1)
	 */
	public double getUniqueKMerCountMatchProb(StatsPerTaxid stats) {
		double mean = getExpectedUniqueKMers(stats);
		double sd = Math.sqrt(getUniqueKMersVariance(stats));
		
		NormalDistribution nd = new NormalDistribution(mean, sd);
		
		return nd.cumulativeProbability(stats.getUniqueKMers()) * 2;
	}

	/*
	 * From https://arxiv.org/pdf/1602.05822.pdf (9)
	 */
	public boolean isProbEstimateInRange(StatsPerTaxid stats) {
		double N = taxidToStoreStats.get(stats.getTaxid()).assignedKMers;
		double A = stats.getKMers();

		if (N <= 5 || A <= 5) {
			return false;
		}

		return 1.4 * Math.pow(N, 0.67) <= A && A <= 1.13 * Math.pow(N, 1.19);
	}
}
