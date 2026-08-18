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
import java.nio.charset.StandardCharsets;

import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.finertree.FTGoalKey;
import org.metagene.genestrip.finertree.FTProject;
import org.metagene.genestrip.make.FileGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.store.Database;
import org.metagene.genestrip.tax.SmallTaxTree;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.*;

/**
 * Writes a CSV file with the per-tax-id database-quality metrics (tp, tp+fp, tp+fn, precision, recall
 * and their weighted averages) computed by {@link DBQualityCountsGoal}, ordered by the taxonomy tree.
 *
 * @param <P> the concrete FT project type
 */
public class DBQualityCSVGoal<P extends FTProject> extends FileGoal<P> {
    private static final DecimalFormat DF = new DecimalFormat("0.00000000", new DecimalFormatSymbols(Locale.US));

    private final ObjectGoal<Database, P> storeGoal;
    private final ObjectGoal<Map<String, DBQualityCountsGoal.Counts>, P> kmersPerTaxGoal;

    /**
     * Creates the goal writing the per-tax-id database-quality CSV file.
     *
     * @param project         the FT project
     * @param key             the goal key
     * @param storeGoal       the goal providing the database (for its taxonomy tree)
     * @param kmersPerTaxGoal the goal providing the per-tax-id quality counts
     * @param deps            further goals this goal depends on
     */
    public DBQualityCSVGoal(P project, FTGoalKey key, ObjectGoal<Database, P> storeGoal, ObjectGoal<Map<String, DBQualityCountsGoal.Counts>, P> kmersPerTaxGoal, Goal<P>... deps) {
        super(project, key, Goal.append(deps, storeGoal, kmersPerTaxGoal));
        this.storeGoal = storeGoal;
        this.kmersPerTaxGoal = kmersPerTaxGoal;
    }

    /**
     * Returns the single CSV output file written by this goal.
     *
     * @return the list containing the CSV output file
     */
    @Override
    public List<File> getFiles() {
        return Collections.singletonList(getProject().getOutputFile(getKey().getName(), GSProject.GSFileType.CSV, false));
    }

    /**
     * Writes the per-tax-id database-quality metrics to the given CSV file, ordered by the taxonomy
     * tree.
     *
     * @param file the CSV file to write
     * @throws IOException if writing the file fails
     */
    @Override
    protected void makeFile(File file) throws IOException {
        Map<String, DBQualityCountsGoal.Counts> kmersPerTax = kmersPerTaxGoal.get();
        SmallTaxTree tree = storeGoal.get().getTaxTree();

        try (PrintStream ps = new PrintStream(file, StandardCharsets.UTF_8)) {
            // The tp/tp+fp/tp+fn columns refer to the paths from the data taxa under a node up to the
            // root and back the unweighted averages and the recalls. The subtree columns refer to the
            // paths up to the node itself and back the weighted avg precision, which is relative to the
            // subtree it is reported for and is simply their quotient. The subtree precision refers to
            // the same subtree and averages the per-k-mer precision over its "subtree kmers" k-mers.
            // The last two columns restrict the subtree precision to the k-mers stored above the data
            // taxa. A k-mer at a data taxon is fixed at a precision of one by construction, and such
            // k-mers are the bulk of a database, so they dominate the unrestricted average while being
            // incapable of improvement. They are appended rather than inserted so that the position of
            // the existing columns stays as it was.
            ps.println("taxid;name;rank;parent taxid;tp;tp+fp;tp+fn;subtree tp;subtree tp+fp;subtree kmers;unweighted avg precision;unweighted avg recall;weighted avg precision;weighted avg recall;node precision;subtree precision;subtree kmers above data;restricted subtree precision");
            // We want result in order of the tree:
            for (SmallTaxTree.SmallTaxIdNode node : tree) {
                DBQualityCountsGoal.Counts counts = kmersPerTax.get(node.getTaxId());
                // Nodes without any data taxon underneath - in particular the artificial "OTHER" nodes
                // introduced by the refinement - carry no genomic evidence at all. Precision and recall
                // are undefined for them, so they are left out entirely instead of reporting NaN.
                if (counts != null && counts.getLeaves() > 0) {
                    ps.print(node.getTaxId());
                    ps.print(";");
                    ps.print(node.getName());
                    ps.print(";");
                    ps.print(node.getRank() == null ? "null" : node.getRank().getName());
                    ps.print(";");
                    SmallTaxTree.SmallTaxIdNode parent = node.getParent();
                    ps.print(parent == null ? "null" : parent.getTaxId());
                    ps.print(";");
                    ps.print(counts.getTp());
                    ps.print(";");
                    ps.print(counts.getTpPlusFp());
                    ps.print(";");
                    ps.print(counts.getTpPlusFn());
                    ps.print(";");
                    ps.print(counts.getSubtreeTp());
                    ps.print(";");
                    ps.print(counts.getSubtreeTpPlusFp());
                    ps.print(";");
                    ps.print(counts.getSubtreeKmerSum());
                    ps.print(";");
                    ps.print(format(counts.getAvgPrecision()));
                    ps.print(";");
                    ps.print(format(counts.getAvgRecall()));
                    ps.print(";");
                    ps.print(format(counts.getPrecision()));
                    ps.print(";");
                    ps.print(format(counts.getRecall()));
                    ps.print(";");
                    ps.print(format(counts.getNodePrecision()));
                    ps.print(";");
                    ps.print(format(counts.getSubtreePrecision()));
                    ps.print(";");
                    ps.print(counts.getSubtreeKMersAboveData());
                    ps.print(";");
                    ps.print(format(counts.getRestrictedSubtreePrecision()));
                    ps.print(";");
                    ps.println();
                }
            }
        }
    }

    /**
     * Formats a metric for the CSV, writing an empty field for an undefined value. Node precision is
     * undefined for nodes without k-mers, and the subtree averages are undefined for subtrees without
     * any, which {@link DBQualityCountsGoal.Counts} reports as {@code NaN}.
     *
     * @param value the metric to format
     * @return the formatted value, or the empty string if it is undefined
     */
    private static String format(double value) {
        return Double.isNaN(value) ? "" : DF.format(value);
    }
}
