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
package org.metagene.genestrip.goals.refseq;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Locale;
import java.util.Map;

import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.make.FileGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.GoalKey;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.tax.Rank;

/**
 * Writes a CSV file summarizing the per-species k-mer rank statistics computed by
 * {@link KMerRankStatsGoal}. The taxonomic ranks are grouped into four intervals and, for each
 * interval, the mean number of a species' k-mers falling into that interval (averaged over all
 * species) is written together with the (population) standard deviation over all species, as well as
 * the box-plot quartiles (first quartile {@code q1}, {@code median} and third quartile {@code q3}) of
 * the per-species interval sums.
 * <p>
 * The rank intervals are:
 * <ul>
 * <li>{@code above phylum} - all ranks strictly above {@link Rank#PHYLUM};</li>
 * <li>{@code phylum to genus} - from {@link Rank#PHYLUM} (inclusive) up to {@link Rank#GENUS}
 * (exclusive);</li>
 * <li>{@code genus to species} - from {@link Rank#GENUS} (inclusive) up to {@link Rank#SPECIES}
 * (exclusive);</li>
 * <li>{@code species and below} - {@link Rank#SPECIES} and everything below it. {@link Rank#SPECIES_GROUP}
 * is counted here as well, even though it sits between genus and species by level.</li>
 * </ul>
 * Ranks without a well-defined level ({@link Rank#isIndeterminate()}, e.g. {@code clade} and
 * {@code no rank}) cannot be placed on the phylum/genus/species axis by their rank alone; k-mers
 * landing on such nodes are already resolved to a well-defined representative rank (by tree position)
 * upstream in {@link KMerRankStatsGoal}, so the count arrays reaching this goal only hold well-defined
 * ranks. The indeterminate case below is therefore only a defensive fallback.
 *
 * @param <P> the project type
 */
public class KMerRankStatsCSVGoal<P extends GSProject> extends FileGoal<P> {
	private static final DecimalFormat DF = new DecimalFormat("0.00000000", new DecimalFormatSymbols(Locale.US));

	/** Display names of the rank intervals, indexed as the interval indices below. */
	private static final String[] INTERVAL_NAMES = { "above phylum", "phylum to genus", "genus to species",
			"species and below" };
	private static final int INTERVAL_COUNT = INTERVAL_NAMES.length;

	// Maps each Rank ordinal (the index used in KMerRankStatsGoal's count arrays) to its interval index,
	// or -1 for ranks that have no place on the phylum/genus/species axis (indeterminate ranks).
	private static final int[] INTERVAL_BY_ORDINAL = computeIntervalByOrdinal();

	private final ObjectGoal<Map<String, long[]>, P> rankStatsGoal;

	/**
	 * Creates the goal writing the rank-interval summary CSV file.
	 *
	 * @param project       the project this goal belongs to
	 * @param key           the key identifying this goal
	 * @param rankStatsGoal the goal providing the per-species k-mer rank statistics
	 * @param deps          any further goals this goal depends on
	 */
	@SafeVarargs
	public KMerRankStatsCSVGoal(P project, GoalKey key, ObjectGoal<Map<String, long[]>, P> rankStatsGoal,
								Goal<P>... deps) {
		super(project, key, Goal.append(deps, rankStatsGoal));
		this.rankStatsGoal = rankStatsGoal;
	}

	@Override
	public List<File> getFiles() {
		return Collections.singletonList(
				getProject().getOutputFile(getKey().getName(), GSProject.GSFileType.CSV, false));
	}

	@Override
	protected void makeFile(File file) throws IOException {
		Map<String, long[]> stats = rankStatsGoal.get();

		// Reduce each species' per-rank counts to its four per-interval sums.
		List<long[]> perSpecies = new ArrayList<>(stats.size());
		double[] sum = new double[INTERVAL_COUNT];
		for (long[] counts : stats.values()) {
			long[] intervalSums = new long[INTERVAL_COUNT];
			for (int ordinal = 0; ordinal < counts.length; ordinal++) {
				int interval = INTERVAL_BY_ORDINAL[ordinal];
				if (interval >= 0) {
					intervalSums[interval] += counts[ordinal];
				}
			}
			perSpecies.add(intervalSums);
			for (int i = 0; i < INTERVAL_COUNT; i++) {
				sum[i] += intervalSums[i];
			}
		}

		int n = perSpecies.size();
		double[] mean = new double[INTERVAL_COUNT];
		for (int i = 0; i < INTERVAL_COUNT; i++) {
			mean[i] = n == 0 ? 0 : sum[i] / n;
		}
		// Second pass for the (population) standard deviation - more accurate than a sum-of-squares.
		double[] sqDiff = new double[INTERVAL_COUNT];
		for (long[] intervalSums : perSpecies) {
			for (int i = 0; i < INTERVAL_COUNT; i++) {
				double d = intervalSums[i] - mean[i];
				sqDiff[i] += d * d;
			}
		}
		double[] std = new double[INTERVAL_COUNT];
		for (int i = 0; i < INTERVAL_COUNT; i++) {
			std[i] = n == 0 ? 0 : Math.sqrt(sqDiff[i] / n);
		}

		// The box-plot quartiles (Q1, median, Q3) of the per-species interval sums. For each interval the
		// species values are gathered into their own array and sorted, then the quartiles are read off by
		// linear interpolation (the type-7 / matplotlib / R / NumPy default). These give the box bounds
		// (Q1, Q3) and the median line; the 1.5*IQR whiskers can be derived downstream from Q1 and Q3.
		double[] q1 = new double[INTERVAL_COUNT];
		double[] median = new double[INTERVAL_COUNT];
		double[] q3 = new double[INTERVAL_COUNT];
		for (int i = 0; i < INTERVAL_COUNT; i++) {
			long[] values = new long[n];
			for (int s = 0; s < n; s++) {
				values[s] = perSpecies.get(s)[i];
			}
			Arrays.sort(values);
			q1[i] = quantile(values, 0.25);
			median[i] = quantile(values, 0.5);
			q3[i] = quantile(values, 0.75);
		}

		try (PrintStream ps = new PrintStream(file, StandardCharsets.UTF_8)) {
			ps.println("rank interval;avg kmers per species;stddev;q1;median;q3;species count;");
			for (int i = 0; i < INTERVAL_COUNT; i++) {
				ps.print(INTERVAL_NAMES[i]);
				ps.print(";");
				ps.print(DF.format(mean[i]));
				ps.print(";");
				ps.print(DF.format(std[i]));
				ps.print(";");
				ps.print(DF.format(q1[i]));
				ps.print(";");
				ps.print(DF.format(median[i]));
				ps.print(";");
				ps.print(DF.format(q3[i]));
				ps.print(";");
				ps.print(n);
				ps.print(";");
				ps.println();
			}
		}
	}

	/**
	 * Computes the {@code p}-quantile of an already ascending-sorted array using linear interpolation
	 * between the two closest ranks (the type-7 definition used by NumPy, R and matplotlib's box plots).
	 *
	 * @param sorted the values in ascending order; must not be empty
	 * @param p      the quantile in {@code [0, 1]}, e.g. {@code 0.25} for the first quartile
	 * @return the interpolated quantile value, or {@code 0} if the array is empty
	 */
	protected static double quantile(long[] sorted, double p) {
		int n = sorted.length;
		if (n == 0) {
			return 0;
		}
		if (n == 1) {
			return sorted[0];
		}
		double pos = p * (n - 1);
		int lower = (int) Math.floor(pos);
		double frac = pos - lower;
		if (lower + 1 >= n) {
			return sorted[n - 1];
		}
		return sorted[lower] + frac * (sorted[lower + 1] - sorted[lower]);
	}

	// Precomputes the interval index for every Rank ordinal (see INTERVAL_BY_ORDINAL).
	private static int[] computeIntervalByOrdinal() {
		Rank[] ranks = Rank.values();
		int[] result = new int[ranks.length];
		for (int i = 0; i < ranks.length; i++) {
			result[i] = intervalOf(ranks[i]);
		}
		return result;
	}

	/**
	 * Classifies a rank into one of the four output intervals.
	 *
	 * @param rank the rank to classify
	 * @return the interval index (0-3), or -1 if the rank has no place on the phylum/genus/species axis
	 */
	protected static int intervalOf(Rank rank) {
		if (rank.isIndeterminate()) {
			return -1;
		}
		// The species group sits between genus and species by level but is counted as species here.
		if (rank == Rank.SPECIES_GROUP) {
			return 3;
		}
		int level = rank.getLevel();
		if (level < Rank.PHYLUM.getLevel()) {
			return 0;
		}
		if (level < Rank.GENUS.getLevel()) {
			return 1;
		}
		if (level < Rank.SPECIES.getLevel()) {
			return 2;
		}
		return 3;
	}
}
