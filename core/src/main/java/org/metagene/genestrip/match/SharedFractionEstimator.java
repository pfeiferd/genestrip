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

/**
 * Estimates the share of a tax id's stored k-mers that a sample contains, with a confidence interval,
 * and the evolutionary distance derived from it. The derivations are in {@code SharedFraction.md}.
 * <p>
 * The hits on a stored k-mer the sample shares are Poisson-distributed, so with too few reads some
 * shared k-mers are still unhit and the observed share {@code u / N} understates the true one. Two
 * estimators correct for this:
 * <ul>
 * <li>with one bit per k-mer, {@link #endpoint(long, long, long)} solves
 * {@code u/N = pi (1 - exp(-c / (pi N)))} for {@code pi};</li>
 * <li>with a two-bit counter per k-mer, {@link #chao(long, long, long)} adds Chao's estimate of the
 * unhit shared k-mers, {@code f1 (f1 - 1) / (2 (f2 + 1))}, to {@code u}. It does not use {@code c} and
 * is therefore not thrown off by k-mers hit very often.</li>
 * </ul>
 * The hits of a read fall on neighbouring k-mers, which makes them dependent. The standard errors
 * account for this through {@code b}, the mean length of the run of hits a hit belongs to.
 */
public class SharedFractionEstimator {
	private final int k;
	private final double z;

	/**
	 * Creates an estimator.
	 *
	 * @param k          the k-mer length, for the evolutionary distance
	 * @param confidence the level of the confidence intervals, in {@code (0, 1)}
	 */
	public SharedFractionEstimator(int k, double confidence) {
		if (!(confidence > 0 && confidence < 1)) {
			throw new IllegalArgumentException("confidence must be in (0, 1), got " + confidence);
		}
		this.k = k;
		this.z = normalQuantile(1 - (1 - confidence) / 2);
	}

	/**
	 * Estimates the shared fraction of the given counts and stores it, with its confidence interval and
	 * the evolutionary distance, in the counts. Chao's estimator is used where the k-mers hit once and
	 * twice were counted, the endpoint estimator otherwise. Where the endpoint estimator has no k-mer hit
	 * twice to go by, only the interval {@code [u/n, 1]} is stored and the estimate stays NaN. Counts
	 * without hits or without stored k-mers are left alone.
	 *
	 * @param counts the counts of one tax id, completed (see {@code CountsPerTaxid#completeValues})
	 */
	public void estimate(CountsPerTaxid counts) {
		long n = counts.getDbKMers();
		long c = counts.kmers;
		long u = counts.uniqueKmers;
		if (n <= 0 || c <= 0 || u <= 0) {
			return;
		}
		double b = runLength(counts.contigLenSquaredSum, counts.kmers);
		double[] estimate;
		if (counts.kmersHitOnce >= 0) {
			estimate = chaoInterval(u, counts.kmersHitOnce, counts.kmersHitTwice, n, b);
		} else if (u >= c) {
			// No k-mer was hit twice: the endpoint equation has no solution, and all that is known is
			// that pi lies between u / n and one. Report that interval, but no estimate.
			estimate = new double[] { Double.NaN, Math.min(1, (double) u / n), 1 };
		} else {
			estimate = endpointInterval(c, u, n, b);
		}
		counts.sharedFraction = estimate[0];
		counts.sharedFractionLow = estimate[1];
		counts.sharedFractionHigh = estimate[2];
		counts.evoDistance = distance(estimate[0]);
		counts.evoDistanceLow = distance(estimate[2]);
		counts.evoDistanceHigh = distance(estimate[1]);
	}

	/**
	 * The evolutionary distance {@code 1 - pi^(1/k)}: the share of positions at which two sequences
	 * differ if the missing k-mers were lost to point mutations alone.
	 *
	 * @param pi the shared fraction
	 * @return the evolutionary distance
	 */
	public double distance(double pi) {
		return 1 - Math.pow(pi, 1.0 / k);
	}

	/**
	 * Returns the endpoint estimate with its confidence interval.
	 *
	 * @param c the hits
	 * @param u the k-mers hit at least once
	 * @param n the stored k-mers
	 * @param b the mean run length of hits, see {@link #runLength(long, long)}
	 * @return the estimate, the lower and the upper bound
	 */
	public double[] endpointInterval(long c, long u, long n, double b) {
		double pi = endpoint(c, u, n);
		return interval(pi, endpointStandardError(c, u, n, b), u, n);
	}

	/**
	 * Returns Chao's estimate with its confidence interval.
	 *
	 * @param u  the k-mers hit at least once
	 * @param f1 the k-mers hit exactly once
	 * @param f2 the k-mers hit exactly twice
	 * @param n  the stored k-mers
	 * @param b  the mean run length of hits, see {@link #runLength(long, long)}
	 * @return the estimate, the lower and the upper bound
	 */
	public double[] chaoInterval(long u, long f1, long f2, long n, double b) {
		double pi = Math.min(1, chao(u, f1, f2) / n);
		return interval(pi, chaoStandardError(u, f1, f2, n, b), u, n);
	}

	private double[] interval(double pi, double se, long u, long n) {
		double floor = Math.min(1, (double) u / n);
		return new double[] { pi, Math.max(floor, pi - z * se), Math.min(1, pi + z * se) };
	}

	/**
	 * The mean length of the run of hits a hit belongs to, {@code E[B^2] / E[B]}, from Genestrip's
	 * per-tax-id contig counters: the sum of the squared contig lengths and the sum of the contig
	 * lengths, which is the number of hits. Runs are counted per contig rather than per read, so a read
	 * whose hits are split counts as two runs, which makes {@code b} a little too small.
	 *
	 * @param squaredSum the sum of the squared contig lengths
	 * @param sum        the sum of the contig lengths
	 * @return the mean run length, at least one
	 */
	public static double runLength(long squaredSum, long sum) {
		return sum > 0 ? Math.max(1, (double) squaredSum / sum) : 1;
	}

	/**
	 * Solves {@code u/n = pi (1 - exp(-c / (pi n)))} for {@code pi}, capped at one. The right-hand side
	 * increases in {@code pi}, so there is exactly one solution once some k-mer was hit twice
	 * ({@code u < c}); before that, the estimate is one.
	 *
	 * @param c the hits
	 * @param u the k-mers hit at least once
	 * @param n the stored k-mers
	 * @return the estimated shared fraction
	 */
	public static double endpoint(long c, long u, long n) {
		double y = (double) u / n;
		double x = (double) c / n;
		if (u >= c || y >= 1) {
			return 1;
		}
		if (f(1, x, y) <= 0) {
			return 1;
		}
		double lo = y;
		double hi = 1;
		for (int i = 0; i < 100; i++) {
			double mid = (lo + hi) / 2;
			if (f(mid, x, y) > 0) {
				hi = mid;
			} else {
				lo = mid;
			}
		}
		return (lo + hi) / 2;
	}

	private static double f(double pi, double x, double y) {
		return pi * (1 - Math.exp(-x / pi)) - y;
	}

	/**
	 * The standard error of {@link #endpoint(long, long, long)}: the variance of {@code u} given
	 * {@code c} with runs of {@code b} hits,
	 * {@code M e^-l (1 - (1 + l b) e^-l + 2 sum_{j=1}^{b-1} (e^(l (j/b - 1)) - e^-l))} with
	 * {@code l = c / M}, turned into a standard error of {@code pi} by the delta method.
	 *
	 * @param c the hits
	 * @param u the k-mers hit at least once
	 * @param n the stored k-mers
	 * @param b the mean run length of hits
	 * @return the standard error, infinite where it cannot be determined
	 */
	public static double endpointStandardError(long c, long u, long n, double b) {
		double pi = endpoint(c, u, n);
		double m = pi * n;
		double l = c / m;
		double e = Math.exp(-l);
		double deriv = 1 - (1 + l) * e;
		if (u >= c || deriv <= 0) {
			return Double.POSITIVE_INFINITY;
		}
		int bi = (int) Math.max(1, Math.round(b));
		double sum = 0;
		for (int j = 1; j < bi; j++) {
			// e^-l (e^(l j/b) - 1), written so that it cannot overflow for a large l.
			sum += Math.exp(l * ((double) j / bi - 1)) - e;
		}
		double var = m * e * (1 - (1 + l * bi) * e + 2 * sum);
		return Math.sqrt(Math.max(var, 0)) / (n * deriv);
	}

	/**
	 * Chao's estimate of the number of shared k-mers in its bias-corrected form,
	 * {@code u + f1 (f1 - 1) / (2 (f2 + 1))}, which stays finite when {@code f2 = 0}.
	 *
	 * @param u  the k-mers hit at least once
	 * @param f1 the k-mers hit exactly once
	 * @param f2 the k-mers hit exactly twice
	 * @return the estimated number of shared k-mers
	 */
	public static double chao(long u, long f1, long f2) {
		return u + (double) f1 * (f1 - 1) / (2.0 * (f2 + 1));
	}

	/**
	 * The standard error of Chao's estimate divided by {@code n}: the variance of the estimated unhit
	 * k-mers with runs of {@code b} hits, {@code (f1 / l) sum_{j=-(b-1)}^{b-1} phi(|j|/b)} with
	 * {@code phi(s) = e^(-l s) (1-s)^2 ((1-s)^2 + 4 (1-s)/l + 2/l^2)} and {@code l = 2 f2 / f1}. The
	 * factor {@code f1 / l} stands for {@code M e^-l}, which it equals in expectation; taking it from the
	 * observed {@code f1} keeps the variance small where few k-mers are left unhit, and makes it equal to
	 * Chao's variance for independent draws, {@code f2 (r^4/4 + r^3 + r^2/2)} with {@code r = f1/f2}, for
	 * {@code b = 1}.
	 *
	 * @param u  the k-mers hit at least once
	 * @param f1 the k-mers hit exactly once
	 * @param f2 the k-mers hit exactly twice
	 * @param n  the stored k-mers
	 * @param b  the mean run length of hits
	 * @return the standard error, infinite where nothing was hit twice yet
	 */
	public static double chaoStandardError(long u, long f1, long f2, long n, double b) {
		if (f1 == 0) {
			// No k-mer was hit only once: none is likely to be left unhit.
			return 0;
		}
		if (f2 == 0) {
			return Double.POSITIVE_INFINITY;
		}
		double l = 2.0 * f2 / f1;
		int bi = (int) Math.max(1, Math.round(b));
		double var = phi(0, l);
		for (int j = 1; j < bi; j++) {
			var += 2 * phi((double) j / bi, l);
		}
		return Math.sqrt(f1 / l * var) / n;
	}

	private static double phi(double s, double l) {
		double r = 1 - s;
		return Math.exp(-l * s) * r * r * (r * r + 4 * r / l + 2 / (l * l));
	}

	/**
	 * The quantile function of the standard normal distribution, by Acklam's rational approximation
	 * (relative error below 1.15e-9).
	 *
	 * @param p the probability, in {@code (0, 1)}
	 * @return the quantile
	 */
	public static double normalQuantile(double p) {
		final double[] a = { -3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02,
				1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00 };
		final double[] bb = { -5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02,
				6.680131188771972e+01, -1.328068155288572e+01 };
		final double[] cc = { -7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00,
				-2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00 };
		final double[] d = { 7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00,
				3.754408661907416e+00 };
		final double low = 0.02425;
		if (p < low) {
			double q = Math.sqrt(-2 * Math.log(p));
			return (((((cc[0] * q + cc[1]) * q + cc[2]) * q + cc[3]) * q + cc[4]) * q + cc[5])
					/ ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1);
		}
		if (p > 1 - low) {
			return -normalQuantile(1 - p);
		}
		double q = p - 0.5;
		double r = q * q;
		return (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q
				/ (((((bb[0] * r + bb[1]) * r + bb[2]) * r + bb[3]) * r + bb[4]) * r + 1);
	}
}
