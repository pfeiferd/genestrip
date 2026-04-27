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

import org.metagene.genestrip.finertree.FTConfigKey;
import org.metagene.genestrip.finertree.FTGoalKey;
import org.metagene.genestrip.finertree.FTProject;
import org.metagene.genestrip.finertree.cluster.DendrogramNode;
import org.metagene.genestrip.finertree.cluster.Similarity;
import org.metagene.genestrip.finertree.cluster.SimpleAggloClustering;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.tax.SmallTaxTree;

import java.util.HashMap;
import java.util.Map;

public class DendrogramGoal<P extends FTProject> extends ObjectGoal<Map<SmallTaxTree.SmallTaxIdNode, DendrogramNode>, P> {
    private final ObjectGoal<KMerIntersectCountGoal.IntersectionsPerNode, P> kmerIntersectGoal;

    @SafeVarargs
    public DendrogramGoal(P project, ObjectGoal<KMerIntersectCountGoal.IntersectionsPerNode, P> kmerIntersectGoal, Goal<P>... deps) {
        super(project, FTGoalKey.DENDROGRAM, append(deps, kmerIntersectGoal));
        this.kmerIntersectGoal = kmerIntersectGoal;
    }

    @Override
    protected void doMakeThis() {
        Map<SmallTaxTree.SmallTaxIdNode, DendrogramNode> res = new HashMap<>();
        KMerIntersectCountGoal.IntersectionsPerNode intersections = kmerIntersectGoal.get();
        SimpleAggloClustering.Method method = (SimpleAggloClustering.Method) configValue(FTConfigKey.CLUSTER_METHOD);
        SimpleAggloClustering clustering = new SimpleAggloClustering(method);
        boolean withChildCounts = booleanConfigValue(FTConfigKey.WITH_DESCENDANT_COUNTS);
        for (SmallTaxTree.SmallTaxIdNode parent : intersections.getParentNodes()) {
            DendrogramNode node = clustering.cluster(new Similarity() {
                @Override
                public int values() {
                    int otherPos = parent.getSubNodes().length;
                    long count = intersections.getIntersectionCount(parent, otherPos, otherPos);
                    return  parent.getSubNodes().length + (count == 0 ? 0 : 1); // "+ 1" for "OTHER_VALUE"
                }

                @Override
                public double getSimilarity(int i, int j) {
                    return intersections.getJaccardIndex(parent, i, j, withChildCounts);
                }
            });
            res.put(parent, node);
        }
        set(res);
    }
}
