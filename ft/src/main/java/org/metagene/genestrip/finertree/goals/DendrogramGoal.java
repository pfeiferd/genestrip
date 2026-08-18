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

/**
 * Builds, for every refinement parent node, a {@link DendrogramNode} tree obtained by agglomerative
 * clustering of the parent's child taxa (plus the optional {@code OTHER} pseudo-item) according to
 * their pairwise k-mer Jaccard similarities.
 * <p>
 * The similarities are taken from the {@link KMerIntersectCountGoal.IntersectionsPerNode} result and
 * clustered with {@link SimpleAggloClustering} using the configured
 * {@link FTConfigKey#CLUSTER_METHOD linkage method}. The resulting map associates each
 * refinement-parent {@link SmallTaxTree.SmallTaxIdNode} with the root of its dendrogram.
 *
 * @param <P> the concrete {@link FTProject} type this goal operates on
 */
public class DendrogramGoal<P extends FTProject> extends ObjectGoal<Map<SmallTaxTree.SmallTaxIdNode, DendrogramNode>, P> {
    private final ObjectGoal<KMerIntersectCountGoal.IntersectionsPerNode, P> kmerIntersectGoal;

    /**
     * Creates the dendrogram goal.
     *
     * @param project           the finer-tree project this goal belongs to
     * @param kmerIntersectGoal the goal providing the per-parent k-mer intersection counts used to
     *                          derive pairwise Jaccard similarities
     * @param deps              additional goals this goal depends on
     */
    @SafeVarargs
    public DendrogramGoal(P project, ObjectGoal<KMerIntersectCountGoal.IntersectionsPerNode, P> kmerIntersectGoal, Goal<P>... deps) {
        super(project, FTGoalKey.DENDROGRAM, append(deps, kmerIntersectGoal));
        this.kmerIntersectGoal = kmerIntersectGoal;
    }

    /**
     * Computes the dendrogram for each refinement parent and stores the resulting map as this goal's
     * value. For every parent the child taxa are clustered by their pairwise
     * {@link KMerIntersectCountGoal.IntersectionsPerNode#getJaccardIndex Jaccard indices}; the
     * {@code OTHER} pseudo-item is only included as an additional leaf when it actually has k-mer
     * intersections (a non-zero self-count).
     */
    @Override
    protected void doMakeThis() {
        Map<SmallTaxTree.SmallTaxIdNode, DendrogramNode> res = new HashMap<>();
        KMerIntersectCountGoal.IntersectionsPerNode intersections = kmerIntersectGoal.get();
        SimpleAggloClustering.Method method = (SimpleAggloClustering.Method) configValue(FTConfigKey.CLUSTER_METHOD);
        SimpleAggloClustering clustering = new SimpleAggloClustering(method);
        boolean withDescendantCounts = booleanConfigValue(FTConfigKey.WITH_DESCENDANT_COUNTS);
        for (SmallTaxTree.SmallTaxIdNode parent : intersections.getParentNodes()) {
            DendrogramNode node = clustering.cluster(new Similarity() {
                @Override
                public int values() {
                    int otherPos = parent.getSubNodes().length;
                    long count = intersections.getIntersectionCount(parent, otherPos, otherPos, withDescendantCounts);
                    return  parent.getSubNodes().length + (count == 0 ? 0 : 1); // "+ 1" for "OTHER_VALUE"
                }

                @Override
                public double getSimilarity(int i, int j) {
                    return intersections.getJaccardIndex(parent, i, j, withDescendantCounts);
                }
            });
            res.put(parent, node);
        }
        set(res);
    }
}
