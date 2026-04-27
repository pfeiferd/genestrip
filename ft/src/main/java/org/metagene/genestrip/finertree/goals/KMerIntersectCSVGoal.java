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
import org.metagene.genestrip.GSConfigKey;
import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.finertree.FTConfigKey;
import org.metagene.genestrip.finertree.FTGoalKey;
import org.metagene.genestrip.finertree.FTProject;
import org.metagene.genestrip.io.StreamProvider;
import org.metagene.genestrip.make.*;
import org.metagene.genestrip.store.Database;
import org.metagene.genestrip.tax.SmallTaxTree;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.*;

public class KMerIntersectCSVGoal<P extends FTProject> extends FileListGoal<P> {
    private static final DecimalFormat DF = new DecimalFormat("0.00000000", new DecimalFormatSymbols(Locale.US));

    private final ObjectGoal<Database, P> storeGoal;
    private final ObjectGoal<KMerIntersectCountGoal.IntersectionsPerNode, P> kmerIntersectGoal;
    private final Map<File, SmallTaxTree.SmallTaxIdNode> fileToNodeMap;

    @SafeVarargs
    public KMerIntersectCSVGoal(P project, ObjectGoal<Database, P> storeGoal, ObjectGoal<KMerIntersectCountGoal.IntersectionsPerNode, P> kmerIntersectGoal, Goal<P>... deps) {
        super(project, FTGoalKey.INTERSECT_CSV, (List<File>) null,  append(deps, kmerIntersectGoal));
        this.storeGoal = storeGoal;
        this.kmerIntersectGoal = kmerIntersectGoal;
        fileToNodeMap = new HashMap<>();
    }

    @Override
    public boolean isCleaned() {
        // Actual check would require make of other goals - not doable.
        return true;
    }

    @Override
    // Do not access kmerIntersectGoal here as it would trigger the related computation already...
    protected void provideFiles() {
        Database database = storeGoal.get();
        SmallTaxTree tree = database.getTaxTree();
        Collection<SmallTaxTree.SmallTaxIdNode> parents = getNodesForPositions(database, (Collection<FTConfigKey.RefinementPosition>) configValue(FTConfigKey.REFINEMENT_POSITIONS));
        List<String> taxids = (List<String>) configValue(GSConfigKey.TAX_IDS);
        for (String taxid : taxids) {
            SmallTaxTree.SmallTaxIdNode parentNode = tree.getNodeByTaxId(taxid);
            if (parentNode != null) {
                parents.add(parentNode);
            }
        }
        for (SmallTaxTree.SmallTaxIdNode node : parents) {
            File matchFile = getProject().getOutputFile(getKey().getName(), node.getTaxId(), null, GSProject.GSFileType.CSV, false);
            addFile(matchFile);
            fileToNodeMap.put(matchFile, node);
        }
    }

    @Override
    protected void makeFile(File file) throws IOException {
        SmallTaxTree.SmallTaxIdNode parent = fileToNodeMap.get(file);
        KMerIntersectCountGoal.IntersectionsPerNode intersections = kmerIntersectGoal.get();

        try (PrintStream out = new PrintStream(StreamProvider.getOutputStreamForFile(file))) {
            out.println("children; kmer sum; kmer spread sum; avg kmer spread; overspread ratio;");
            int nChildren = parent.getSubNodes().length;
            out.print(nChildren);
            out.print(';');
            out.print(intersections.getKMerSum(parent));
            out.print(';');
            out.print(intersections.getKMerSpreadSum(parent));
            out.print(';');
            out.print(DF.format(intersections.getAvgKMerSpread(parent)));
            out.print(';');
            out.print(DF.format(intersections.getOverspreadRatio(parent)));
            out.print(';');
            out.println();

            SmallTaxTree.SmallTaxIdNode[] children = parent.getSubNodes();
            for (int i = 0; i < children.length; i++) {
                out.print(children[i].getTaxId());
                out.print(';');
            }
            out.print("OTHER;");
            out.println();
            for (int i = 0; i <= children.length; i++) {
                for (int j = 0; j <= children.length; j++) {
                    out.print(intersections.getIntersectionCount(parent, i, j));
                    out.print(';');
                }
                out.println();
            }
            out.println();
            for (int i = 0; i <= children.length; i++) {
                for (int j = 0; j <= children.length; j++) {
                    out.print(DF.format(intersections.getJaccardIndex(parent, i, j, booleanConfigValue(FTConfigKey.WITH_DESCENDANT_COUNTS))));
                    out.print(';');
                }
                out.println();
            }
        }
    }

    public static Collection<SmallTaxTree.SmallTaxIdNode> getNodesForPositions(Database database, Collection<FTConfigKey.RefinementPosition> refinementIntervals) {
        Collection<SmallTaxTree.SmallTaxIdNode> res = new ArrayList<>();
        Object2LongMap<String> stats = database.getStats();
        Iterator<SmallTaxTree.SmallTaxIdNode> it = database.getTaxTree().iterator();
        while (it.hasNext()) {
            SmallTaxTree.SmallTaxIdNode node = it.next();
            if (FTConfigKey.RefinementPosition.getMatchingNodeFor(node, refinementIntervals) != null) {
                if (node.getSubNodes() != null && node.getSubNodes().length > 0) {
                    if (stats.getOrDefault(node.getTaxId(), 0) > 0) {
                        res.add(node);
                    }
                }
            }
        }
        return res;
    }
}
