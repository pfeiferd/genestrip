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
package org.metagene.genestrip.finertree.goals;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;

import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.finertree.FTGoalKey;
import org.metagene.genestrip.finertree.FTProject;
import org.metagene.genestrip.make.FileGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.store.Database;
import org.metagene.genestrip.tax.Rank;
import org.metagene.genestrip.tax.SmallTaxTree;
import org.metagene.genestrip.tax.SmallTaxTree.SmallTaxIdNode;

/**
 * Writes a single CSV file that aggregates the per-node branching-degree histograms from
 * {@link KMerBranchHistoGoal} by taxonomic rank, with one row per rank.
 * <p>
 * Every node is first turned into a relative branching-degree distribution: its k-mer count per
 * branching degree divided by the node's total k-mer count, with the trailing "OTHER" column (the last
 * histogram slot) excluded from both the total and the degrees, so the relative values of a node sum to
 * one. Nodes are then grouped by rank and, for each branching degree from {@code 1} up to
 * {@link #MAX_DEGREE}, the relative values of all nodes of that rank that actually have that degree
 * (i.e. that have at least that many child nodes) form a sample. For each such (rank, branching degree)
 * pair the sample's mean, (population) standard deviation and first quartile, median and third quartile
 * (type-7 quantiles, matching NumPy / R / matplotlib) are written, laid out along the row and named by
 * combining the branching degree with the statistic (e.g. {@code 1-avg}, {@code 1-stddev}, ...,
 * {@code 10-q3}). Branching degrees beyond {@link #MAX_DEGREE} are not reported; degrees a rank never
 * reaches are left empty.
 * <p>
 * Each row additionally starts with two summary blocks over a per-node branching-degree measure:
 * {@code childdeg-*} over the node's number of children in the database's taxonomy tree (structural,
 * excluding the OTHER column, independent of the k-mers), and {@code kmerdeg-*} over the k-mer-weighted
 * mean branching degree of the node (the weighted mean of the branching degree over its k-mers, using its
 * full distribution, uncapped, excluding the OTHER column). Since these two blocks report branching
 * degrees themselves rather than relative k-mer shares, they carry the sample's minimum and maximum on
 * top of the five statistics above, written as {@code avg, stddev, min, q1, median, q3, max}, so that the
 * extremes enclose the quartiles.
 * Both measures are summarised across the rank's nodes. Rows are ordered by rank.
 *
 * @param <P> the concrete FT project type
 */
public class KMerBranchHistoRankCSVGoal<P extends FTProject> extends FileGoal<P> {
    /** Highest branching degree for which statistics are written. */
    public static final int MAX_DEGREE = 10;

    private static final DecimalFormat DF = new DecimalFormat("0.00000000", new DecimalFormatSymbols(Locale.US));

    private final ObjectGoal<Database, P> storeGoal;
    private final ObjectGoal<Map<String, long[]>, P> branchHistoGoal;

    /**
     * Creates the goal under the {@link FTGoalKey#BRANCH_HISTO_RANK_CSV} key.
     *
     * @param project         the FT project
     * @param storeGoal       the goal providing the loaded database (used to resolve each node's rank)
     * @param branchHistoGoal the goal providing the per-node branching-degree histograms
     * @param deps            further goals this goal depends on
     */
    @SafeVarargs
    public KMerBranchHistoRankCSVGoal(P project, ObjectGoal<Database, P> storeGoal,
                                      ObjectGoal<Map<String, long[]>, P> branchHistoGoal, Goal<P>... deps) {
        super(project, FTGoalKey.BRANCH_HISTO_RANK_CSV, Goal.append(deps, storeGoal, branchHistoGoal));
        this.storeGoal = storeGoal;
        this.branchHistoGoal = branchHistoGoal;
    }

    @Override
    public List<File> getFiles() {
        return Collections.singletonList(
                getProject().getOutputFile(getKey().getName(), GSProject.GSFileType.CSV, false));
    }

    @Override
    protected void makeFile(File file) throws IOException {
        Map<String, long[]> histos = branchHistoGoal.get();
        SmallTaxTree tree = storeGoal.get().getTaxTree();

        // Per rank, one growable list of samples indexed by (branching degree - 1), capped at
        // MAX_DEGREE. Each sample holds the relative frequency of that branching degree, one value per
        // node of the rank that has that degree (i.e. that has at least that many children).
        Map<Rank, List<List<Double>>> samplesByRank = new HashMap<>();
        // Per rank, two per-node branching-degree measures: childdeg = the node's number of children in
        // the database's taxonomy tree (structural, excluding the OTHER column); kmerdeg = the k-mer-
        // weighted mean branching degree over the node's k-mers (using its full distribution, uncapped,
        // excluding the OTHER column).
        Map<Rank, List<Double>> childDegByRank = new HashMap<>();
        Map<Rank, List<Double>> kmerDegByRank = new HashMap<>();
        for (Map.Entry<String, long[]> e : histos.entrySet()) {
            SmallTaxIdNode node = tree.getNodeByTaxId(e.getKey());
            Rank rank = node == null ? null : node.getRank();
            long[] histo = e.getValue();
            // Degrees 1..L-1 (indices 0..L-2); the last slot (index L-1) is the OTHER column, excluded.
            // L-1 is therefore also the node's number of children in the tree (its branching degree).
            int degrees = histo.length - 1;
            if (degrees <= 0) {
                continue;
            }
            long total = 0;
            double weighted = 0;
            for (int i = 0; i < degrees; i++) {
                total += histo[i];
                weighted += (double) (i + 1) * histo[i];
            }
            if (total == 0) {
                continue;
            }
            childDegByRank.computeIfAbsent(rank, k -> new ArrayList<>()).add((double) degrees);
            kmerDegByRank.computeIfAbsent(rank, k -> new ArrayList<>()).add(weighted / total);
            int cap = Math.min(degrees, MAX_DEGREE);
            List<List<Double>> perDegree = samplesByRank.computeIfAbsent(rank, k -> new ArrayList<>());
            while (perDegree.size() < cap) {
                perDegree.add(new ArrayList<>());
            }
            for (int i = 0; i < cap; i++) {
                perDegree.get(i).add((double) histo[i] / total);
            }
        }

        List<Rank> ranks = new ArrayList<>(samplesByRank.keySet());
        ranks.sort(Comparator.comparingInt(r -> r == null ? Integer.MAX_VALUE : r.ordinal()));

        try (PrintStream ps = new PrintStream(file, StandardCharsets.UTF_8)) {
            StringBuilder header = new StringBuilder("rank;");
            header.append("childdeg-avg;childdeg-stddev;childdeg-min;childdeg-q1;childdeg-median;childdeg-q3;childdeg-max;");
            header.append("kmerdeg-avg;kmerdeg-stddev;kmerdeg-min;kmerdeg-q1;kmerdeg-median;kmerdeg-q3;kmerdeg-max;");
            for (int d = 1; d <= MAX_DEGREE; d++) {
                header.append(d).append("-avg;").append(d).append("-stddev;").append(d).append("-q1;")
                        .append(d).append("-median;").append(d).append("-q3;");
            }
            ps.println(header.toString());

            for (Rank rank : ranks) {
                List<List<Double>> perDegree = samplesByRank.get(rank);
                ps.print(rank == null ? "no rank" : rank.getName());
                ps.print(';');
                // Summary blocks: distribution across the rank's nodes of the number-of-children branching
                // degree (childdeg) and of the k-mer-weighted average branching degree (kmerdeg).
                printStats(ps, childDegByRank.get(rank), true);
                printStats(ps, kmerDegByRank.get(rank), true);
                for (int d = 1; d <= MAX_DEGREE; d++) {
                    // A degree the rank never reaches -> printStats writes five empty cells.
                    List<Double> values = (d - 1) < perDegree.size() ? perDegree.get(d - 1) : null;
                    printStats(ps, values);
                }
                ps.println();
            }
        }
    }

    /**
     * Computes and writes the five statistics (mean, population standard deviation, and the type-7 first
     * quartile, median and third quartile) of the given sample as five {@code ;}-terminated cells. When
     * the sample is empty or {@code null}, five empty cells are written instead.
     *
     * @param ps     the stream to write to
     * @param values the sample values (unsorted), or {@code null}/empty for five empty cells
     */
    private void printStats(PrintStream ps, List<Double> values) {
        printStats(ps, values, false);
    }

    /**
     * Computes and writes the statistics of the given sample as {@code ;}-terminated cells: mean,
     * population standard deviation, and the type-7 first quartile, median and third quartile, optionally
     * bracketed by the sample's minimum and maximum. The cells are written in the order
     * {@code avg, stddev, [min,] q1, median, q3 [, max]}, so that the minimum and maximum enclose the
     * quartiles and thereby complete the five-number summary. When the sample is empty or {@code null},
     * the corresponding number of empty cells is written instead.
     *
     * @param ps         the stream to write to
     * @param values     the sample values (unsorted), or {@code null}/empty for empty cells
     * @param withRange  whether to additionally write the sample's minimum and maximum
     */
    private void printStats(PrintStream ps, List<Double> values, boolean withRange) {
        if (values == null || values.isEmpty()) {
            ps.print(withRange ? ";;;;;;;" : ";;;;;");
            return;
        }
        double[] sample = toSortedArray(values);
        int n = sample.length;
        double mean = 0;
        for (double v : sample) {
            mean += v;
        }
        mean /= n;
        double var = 0;
        for (double v : sample) {
            double diff = v - mean;
            var += diff * diff;
        }
        double std = Math.sqrt(var / n);
        ps.print(DF.format(mean));
        ps.print(';');
        ps.print(DF.format(std));
        ps.print(';');
        if (withRange) {
            // The sample is sorted ascending, so the extremes are its first and last element.
            ps.print(DF.format(sample[0]));
            ps.print(';');
        }
        ps.print(DF.format(quantile(sample, 0.25)));
        ps.print(';');
        ps.print(DF.format(quantile(sample, 0.5)));
        ps.print(';');
        ps.print(DF.format(quantile(sample, 0.75)));
        ps.print(';');
        if (withRange) {
            ps.print(DF.format(sample[n - 1]));
            ps.print(';');
        }
    }

    /**
     * Copies the given values into an ascending-sorted array.
     *
     * @param values the values to sort
     * @return a new array holding the values in ascending order
     */
    private static double[] toSortedArray(List<Double> values) {
        double[] a = new double[values.size()];
        for (int i = 0; i < a.length; i++) {
            a[i] = values.get(i);
        }
        java.util.Arrays.sort(a);
        return a;
    }

    /**
     * Computes the {@code p}-quantile of an already ascending-sorted array using linear interpolation
     * between the two closest ranks (the type-7 definition used by NumPy, R and matplotlib).
     *
     * @param sorted the values in ascending order
     * @param p      the quantile in {@code [0, 1]}
     * @return the interpolated quantile, or {@code 0} if the array is empty
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
        int lo = (int) Math.floor(pos);
        double frac = pos - lo;
        if (lo + 1 >= n) {
            return sorted[n - 1];
        }
        return sorted[lo] + frac * (sorted[lo + 1] - sorted[lo]);
    }
}
