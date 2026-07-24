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

import org.metagene.genestrip.GSConfigKey;
import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.make.FileGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.GoalKey;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.store.Database;
import org.metagene.genestrip.tax.Rank;
import org.metagene.genestrip.tax.SmallTaxTree;
import org.metagene.genestrip.tax.SmallTaxTree.SmallTaxIdNode;

/**
 * Writes a CSV file summarizing the per-species k-mer rank statistics computed by
 * {@link KMerRankStatsGoal}. The taxonomic ranks are grouped into four intervals and, for each
 * interval, the mean number of a species' k-mers falling into that interval (averaged over all
 * species) is written together with the (population) standard deviation over all species, as well as
 * the box-plot quartiles (first quartile {@code q1}, {@code median} and third quartile {@code q3}) of
 * the per-species interval sums. In addition, the same statistics ({@code rel avg}, {@code rel stddev},
 * {@code rel q1}, {@code rel median}, {@code rel q3}) are written for the per-species relative interval
 * values - each species' interval sum divided by its total over the four intervals, so those ratios sum
 * to one per species. The considered set of species is the actual species-rank taxa of the database (as
 * keyed by {@link KMerRankStatsGoal}; higher-rank taxa are not pseudo-species). A species without any
 * k-mers has no defined ratios and does not enter the statistics, so {@code species count} counts the
 * database's species that carry k-mers and applies to both the absolute and the relative columns.
 * <p>
 * By default every species of the database contributes to these figures. If the configuration key
 * {@link GSConfigKey#KMER_RANK_STATS_TAXID} is set to a tax id, only species that are taxonomic
 * descendants of that node (the node itself included) are taken into account.
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
	private final ObjectGoal<Database, P> dbGoal;

	/**
	 * Creates the goal writing the rank-interval summary CSV file.
	 *
	 * @param project       the project this goal belongs to
	 * @param key           the key identifying this goal
	 * @param rankStatsGoal the goal providing the per-species k-mer rank statistics
	 * @param dbGoal        the goal supplying the loaded database, used to resolve the optional tax node
	 *                      the statistics are restricted to (see {@link GSConfigKey#KMER_RANK_STATS_TAXID})
	 * @param deps          any further goals this goal depends on
	 */
	@SafeVarargs
	public KMerRankStatsCSVGoal(P project, GoalKey key, ObjectGoal<Map<String, long[]>, P> rankStatsGoal,
								ObjectGoal<Database, P> dbGoal, Goal<P>... deps) {
		super(project, key, Goal.append(deps, rankStatsGoal, dbGoal));
		this.rankStatsGoal = rankStatsGoal;
		this.dbGoal = dbGoal;
	}

	@Override
	public List<File> getFiles() {
		return Collections.singletonList(
				getProject().getOutputFile(getKey().getName(), GSProject.GSFileType.CSV, false));
	}

	@Override
	protected void makeFile(File file) throws IOException {
		Map<String, long[]> stats = rankStatsGoal.get();

		// Optional restriction: if a tax id is configured, only species that are taxonomic descendants of
		// that node (the node itself included) are included in the statistics. If unset, all species count.
		String taxId = stringConfigValue(GSConfigKey.KMER_RANK_STATS_TAXID);
		SmallTaxTree taxTree = null;
		SmallTaxIdNode filterNode = null;
		if (taxId != null && !taxId.isEmpty()) {
			taxTree = dbGoal.get().getTaxTree();
			filterNode = taxTree.getNodeByTaxId(taxId);
			if (filterNode == null) {
				throw new IllegalStateException("Config key " + GSConfigKey.KMER_RANK_STATS_TAXID.getName()
						+ " refers to tax id " + taxId + " which is not in the database's taxonomy tree.");
			}
		}

		// Reduce each species' per-rank counts to its four per-interval sums. In addition, each species is
		// turned into its relative interval distribution (each interval sum divided by the species' total
		// over the four intervals), so those four ratios sum to one per species. Species with no k-mers at
		// all (total zero) have no defined ratios and are excluded from every statistic here.
		List<long[]> perSpecies = new ArrayList<>(stats.size());
		List<double[]> perSpeciesRel = new ArrayList<>(stats.size());
		double[] sum = new double[INTERVAL_COUNT];
		for (Map.Entry<String, long[]> entry : stats.entrySet()) {
			if (filterNode != null) {
				SmallTaxIdNode speciesNode = taxTree.getNodeByTaxId(entry.getKey());
				if (speciesNode == null || !taxTree.isAncestorOf(speciesNode, filterNode)) {
					continue;
				}
			}
			long[] counts = entry.getValue();
			long[] intervalSums = new long[INTERVAL_COUNT];
			long total = 0;
			for (int ordinal = 0; ordinal < counts.length; ordinal++) {
				int interval = INTERVAL_BY_ORDINAL[ordinal];
				if (interval >= 0) {
					intervalSums[interval] += counts[ordinal];
					total += counts[ordinal];
				}
			}
			if (total == 0) {
				// Species without any k-mers do not enter any statistic (absolute or relative).
				continue;
			}
			perSpecies.add(intervalSums);
			for (int i = 0; i < INTERVAL_COUNT; i++) {
				sum[i] += intervalSums[i];
			}
			double[] rel = new double[INTERVAL_COUNT];
			for (int i = 0; i < INTERVAL_COUNT; i++) {
				rel[i] = (double) intervalSums[i] / total;
			}
			perSpeciesRel.add(rel);
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

		// The same set of statistics over the per-species relative interval values (which sum to one per
		// species). The relative and absolute statistics now share the same species set (those with at
		// least one k-mer), so their per-interval averages sum to one across the four intervals.
		double[] relMean = new double[INTERVAL_COUNT];
		for (double[] rel : perSpeciesRel) {
			for (int i = 0; i < INTERVAL_COUNT; i++) {
				relMean[i] += rel[i];
			}
		}
		for (int i = 0; i < INTERVAL_COUNT; i++) {
			relMean[i] = n == 0 ? 0 : relMean[i] / n;
		}
		double[] relSqDiff = new double[INTERVAL_COUNT];
		for (double[] rel : perSpeciesRel) {
			for (int i = 0; i < INTERVAL_COUNT; i++) {
				double d = rel[i] - relMean[i];
				relSqDiff[i] += d * d;
			}
		}
		double[] relStd = new double[INTERVAL_COUNT];
		for (int i = 0; i < INTERVAL_COUNT; i++) {
			relStd[i] = n == 0 ? 0 : Math.sqrt(relSqDiff[i] / n);
		}
		double[] relQ1 = new double[INTERVAL_COUNT];
		double[] relMedian = new double[INTERVAL_COUNT];
		double[] relQ3 = new double[INTERVAL_COUNT];
		for (int i = 0; i < INTERVAL_COUNT; i++) {
			double[] values = new double[n];
			for (int s = 0; s < n; s++) {
				values[s] = perSpeciesRel.get(s)[i];
			}
			Arrays.sort(values);
			relQ1[i] = quantile(values, 0.25);
			relMedian[i] = quantile(values, 0.5);
			relQ3[i] = quantile(values, 0.75);
		}

		try (PrintStream ps = new PrintStream(file, StandardCharsets.UTF_8)) {
			ps.println("rank interval;avg kmers per species;stddev;q1;median;q3;species count;"
					+ "rel avg;rel stddev;rel q1;rel median;rel q3;");
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
				ps.print(DF.format(relMean[i]));
				ps.print(";");
				ps.print(DF.format(relStd[i]));
				ps.print(";");
				ps.print(DF.format(relQ1[i]));
				ps.print(";");
				ps.print(DF.format(relMedian[i]));
				ps.print(";");
				ps.print(DF.format(relQ3[i]));
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

	/**
	 * Computes the {@code p}-quantile of an already ascending-sorted array using linear interpolation
	 * between the two closest ranks (the type-7 definition used by NumPy, R and matplotlib's box plots).
	 *
	 * @param sorted the values in ascending order
	 * @param p      the quantile in {@code [0, 1]}, e.g. {@code 0.25} for the first quartile
	 * @return the interpolated quantile value, or {@code 0} if the array is empty
	 */
	protected static double quantile(double[] sorted, double p) {
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
