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

public class DBQualityCSVGoal<P extends FTProject> extends FileGoal<P> {
    private static final DecimalFormat DF = new DecimalFormat("0.00000000", new DecimalFormatSymbols(Locale.US));

    private final ObjectGoal<Database, P> storeGoal;
    private final ObjectGoal<Map<String, DBQualityCountsGoal.Counts>, P> kmersPerTaxGoal;

    public DBQualityCSVGoal(P project, FTGoalKey key, ObjectGoal<Database, P> storeGoal, ObjectGoal<Map<String, DBQualityCountsGoal.Counts>, P> kmersPerTaxGoal, Goal<P>... deps) {
        super(project, key, Goal.append(deps, storeGoal, kmersPerTaxGoal));
        this.storeGoal = storeGoal;
        this.kmersPerTaxGoal = kmersPerTaxGoal;
    }

    @Override
    public List<File> getFiles() {
        return Collections.singletonList(getProject().getOutputFile(getKey().getName(), GSProject.GSFileType.CSV, false));
    }

    @Override
    protected void makeFile(File file) throws IOException {
        Map<String, DBQualityCountsGoal.Counts> kmersPerTax = kmersPerTaxGoal.get();
        SmallTaxTree tree = storeGoal.get().getTaxTree();

        try (PrintStream ps = new PrintStream(file, StandardCharsets.UTF_8)) {
            ps.println("taxid;name;rank;parent taxid;tp;tp+fp;tp+fn;precision;recall;weighted avg precision;weighted avg recall;");
            // We want result in order of the tree:
            Iterator<SmallTaxTree.SmallTaxIdNode> iterator = tree.iterator();
            while (iterator.hasNext()) {
                SmallTaxTree.SmallTaxIdNode node = iterator.next();
                DBQualityCountsGoal.Counts counts = kmersPerTax.get(node.getTaxId());
                if (counts != null) {
                    ps.print(node.getTaxId());
                    ps.print(";");
                    ps.print(node.getName());
                    ps.print(";");
                    ps.print(node.getRank().getName());
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
                    ps.print(DF.format(counts.getAvgPrecision()));
                    ps.print(";");
                    ps.print(DF.format(counts.getAvgRecall()));
                    ps.print(";");
                    ps.print(DF.format(counts.getPrecision()));
                    ps.print(";");
                    ps.print(DF.format(counts.getRecall()));
                    ps.print(";");
                    ps.println();
                }
            }
        }
    }
}
