package org.metagene.genestrip.match;

import org.metagene.genestrip.match.FastqKMerMatcher.StatsPerTaxid;
import org.metagene.genestrip.store.KMerStoreWrapper;

import it.unimi.dsi.fastutil.objects.Object2LongMap;

public class UniqueKMerEstimator {
	private long totalKMers;
	private final Object2LongMap<String> storeStats;

	public UniqueKMerEstimator(KMerStoreWrapper storeWrapper) {
		storeStats = storeWrapper.getStats();
	}

	public void setTotalKMers(long totalKMers) {
		this.totalKMers = totalKMers;
	}

	public long getTotalKMers() {
		return totalKMers;
	}

	public double getNormalizedKMers(StatsPerTaxid stats) {
		long totalStoredWithCounts = storeStats.getLong(null);
		return ((double) stats.getKMers()) * totalStoredWithCounts / totalKMers
				/ storeStats.getLong(stats.getTaxid());
	}

	public double getExpectedUniqueKMers(StatsPerTaxid stats) {
		return storeStats.getLong(stats.getTaxid()) * getAtLeastOnceUniqueKMerProb(stats);
	}

	public double getAtLeastOnceUniqueKMerProb(StatsPerTaxid stats) {
		double N = storeStats.getLong(stats.getTaxid());
		double A = stats.getKMers();
		return 1 - Math.pow(1 - 1d / N, A);
	}
	
// Was not useful in practice: Removed for now:	
//	/*
//	 * From https://arxiv.org/pdf/1602.05822.pdf (5)
//	 */
//	public double getUniqueKMersVariance(StatsPerTaxid stats) {
//		double N = storeStats.getLong(stats.getTaxid());
//		double A = stats.getKMers();
//
//		return N * (N - 1) * Math.pow(1 - 2 / N, A) + N * Math.pow(1 - 1 / N, A) - N * N * Math.pow(1 - 1 / N, 2 * A);
//	}
//
//
//	/*
//	 * From https://arxiv.org/pdf/1602.05822.pdf Normal distribution approximates
//	 * P(k) from (1)
//	 */
//	public double getUniqueKMerCountMatchScore(StatsPerTaxid stats) {
//		double mean = getExpectedUniqueKMers(stats);
//		double sd = Math.sqrt(getUniqueKMersVariance(stats));
//
//		if (sd == 0) {
//			return 0;
//		}
//		NormalDistribution nd = new NormalDistribution(mean, sd);
//
//		// We use the densitiy at X = unique kmers. We devide it by the maximum density
//		// (at the mean) to obtain a score value between 0 and 1.
//		return nd.density(stats.getUniqueKMers()) / nd.density(mean);
//
//// Old attempt: P(X <= mean - |mean - unique kmers| or X >= mean + |mean - unique kmers|)  		
////		return nd.cumulativeProbability(mean - Math.abs(mean - stats.getUniqueKMers())) * 2;
//	}
//
//	/*
//	 * From https://arxiv.org/pdf/1602.05822.pdf (9)
//	 */
//	public boolean isProbEstimateInRange(StatsPerTaxid stats) {
//		double N = storeStats.getLong(stats.getTaxid());
//		double A = stats.getKMers();
//
//		if (N <= 5 || A <= 5) {
//			return false;
//		}
//
//		return 1.4 * Math.pow(N, 0.67) <= A && A <= 1.13 * Math.pow(N, 1.19);
//	}
}
