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

import it.unimi.dsi.fastutil.objects.Object2LongMap;
import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.finertree.FTGoalKey;
import org.metagene.genestrip.finertree.FTProject;
import org.metagene.genestrip.make.FileGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.store.Database;
import org.metagene.genestrip.tax.Rank;
import org.metagene.genestrip.tax.SmallTaxTree;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

public class FTConsistencyGoal<P extends FTProject> extends FileGoal<P> {
    private final ObjectGoal<Database, P> storeGoal;
    private final ObjectGoal<Object2LongMap<String>, P> kmersPerTaxGoal;

    public FTConsistencyGoal(P project, FTGoalKey key, ObjectGoal<Database, P> storeGoal, ObjectGoal<Object2LongMap<String>, P> kmersPerTaxGoal, Goal<P>... deps) {
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
        Object2LongMap<String> kmersPerTax = kmersPerTaxGoal.get();
        SmallTaxTree tree = storeGoal.get().getTaxTree();
        Object2LongMap stats = storeGoal.get().getStats();

        try (PrintStream ps = new PrintStream(file)) {
            ps.println("taxid;genome kmers;db kmers;node kmers;");
            Iterator<SmallTaxTree.SmallTaxIdNode> iterator = tree.iterator();
            while (iterator.hasNext()) {
                SmallTaxTree.SmallTaxIdNode node = iterator.next();
                // Leaves where originally stored...
                if (acceptNode(node)) {
                    String key = node.getTaxId();
                    long dbSum = getPathSum(node, stats);
                    long value = kmersPerTax.getOrDefault(key, 0);
                    long nodeKmers = stats.getOrDefault(key, 0);
                    ps.print(key);
                    ps.print(";");
                    ps.print(value);
                    ps.print(";");
                    ps.print(dbSum);
                    ps.print(";");
                    ps.print(nodeKmers);
                    ps.println(";");
                }
            }
        }
    }

    public static boolean acceptNode(SmallTaxTree.SmallTaxIdNode node) {
        return Rank.SPECIES.equals(node.getRank()); // node.getSubNodes() == null || node.getSubNodes().length == 0;
    }

    public static long getPathSum(SmallTaxTree.SmallTaxIdNode node, Object2LongMap stats) {
        int res = 0;
        if (node != null) {
            // We start one above the actual node.
            for (node = node.getParent(); node != null; node = node.getParent()) {
                res += stats.getOrDefault(node.getTaxId(), 0L);
            }
        }
        return res;
    }
}
