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
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;

import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.finertree.FTGoalKey;
import org.metagene.genestrip.finertree.FTProject;
import org.metagene.genestrip.make.FileGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.store.Database;
import org.metagene.genestrip.tax.SmallTaxTree;
import org.metagene.genestrip.tax.SmallTaxTree.SmallTaxIdNode;

/**
 * Writes a single CSV file summarizing the per-node branching-degree histograms computed by
 * {@link KMerBranchHistoGoal}. Each map entry becomes one row: the parent node's tax id and name, its
 * number of child nodes, and then the histogram counts - the count of the node's k-mers with branching
 * degree {@code 1}, with branching degree {@code 2} and so on up to {@code children + 1} (the last
 * degree counting the "OTHER" bucket). Rows are therefore ragged: a node with {@code c} children
 * contributes {@code c + 1} histogram columns. Rows are sorted by tax id for deterministic output.
 *
 * @param <P> the concrete FT project type
 */
public class KMerBranchHistoCSVGoal<P extends FTProject> extends FileGoal<P> {
    private final ObjectGoal<Database, P> storeGoal;
    private final ObjectGoal<Map<String, long[]>, P> branchHistoGoal;

    /**
     * Creates the goal under the {@link FTGoalKey#BRANCH_HISTO_CSV} key.
     *
     * @param project         the FT project
     * @param storeGoal       the goal providing the loaded database (used to resolve node names)
     * @param branchHistoGoal the goal providing the per-node branching-degree histograms
     * @param deps            further goals this goal depends on
     */
    @SafeVarargs
    public KMerBranchHistoCSVGoal(P project, ObjectGoal<Database, P> storeGoal,
                                  ObjectGoal<Map<String, long[]>, P> branchHistoGoal, Goal<P>... deps) {
        super(project, FTGoalKey.BRANCH_HISTO_CSV, Goal.append(deps, storeGoal, branchHistoGoal));
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

        List<String> taxIds = new ArrayList<>(histos.keySet());
        Collections.sort(taxIds);

        try (PrintStream ps = new PrintStream(file, StandardCharsets.UTF_8)) {
            ps.println("tax id;name;children;kmers per branching degree (degree 1 .. children+1, last degree counts the OTHER bucket);");
            for (String taxId : taxIds) {
                long[] histo = histos.get(taxId);
                SmallTaxIdNode node = tree.getNodeByTaxId(taxId);
                ps.print(taxId);
                ps.print(";");
                ps.print(node == null ? "" : node.getName());
                ps.print(";");
                // The histogram has one slot per child plus one for the OTHER bucket.
                ps.print(histo.length - 1);
                ps.print(";");
                for (int i = 0; i < histo.length; i++) {
                    ps.print(histo[i]);
                    ps.print(";");
                }
                ps.println();
            }
        }
    }
}
