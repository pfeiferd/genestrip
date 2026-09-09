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

import java.util.HashMap;
import java.util.Map;

import org.metagene.genestrip.probfilter.ProbFilter;
import org.metagene.genestrip.finertree.FTGoalKey;
import org.metagene.genestrip.finertree.FTProject;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.store.Database;
import org.metagene.genestrip.tax.SmallTaxTree.SmallTaxIdNode;

/**
 * Concrete {@link KMerStoreWorkGoal} that builds, per taxonomy parent node at a refinement rank, a
 * histogram over the <em>branching degree</em> of the node's k-mers as recorded by
 * {@link KMerIndexBloomGoal}. The branching degree of a k-mer is the number of the parent's direct
 * child subtrees it occurs in (as tested against the {@link ProbFilter}). For every
 * visited k-mer the histogram slot for its branching degree is incremented.
 * <p>
 * The trailing "OTHER" bucket - k-mers attributed to no concrete child genome - is treated like an
 * additional child, so the branching degree ranges from {@code 1} up to {@code childCount + 1}.
 * <p>
 * The result is a {@code Map<String, long[]>} keyed by the parent's tax id. In each value array the
 * slot at index {@code branchingDegree - 1} holds the number of the parent's k-mers with that
 * branching degree, so the array is exactly {@code childCount + 1} long (the maximum possible
 * branching degree, counting the OTHER bucket). K-mers that occur in nothing at all (branching
 * degree {@code 0}) are not counted.
 * <p>
 * Only parent nodes that sit at a refinement position and carry k-mers - i.e. the very nodes for
 * which {@link KMerIndexBloomGoal} made entries - end up in the map. Leaf nodes (without children)
 * and nodes without any counted k-mer are absent.
 *
 * @param <P> the concrete FT project type
 */
public class KMerBranchHistoGoal<P extends FTProject> extends KMerStoreWorkGoal<Map<String, long[]>, P> {
    private Map<String, long[]> map;

    /**
     * Creates the goal under the {@link FTGoalKey#BRANCH_HISTO} key.
     *
     * @param project         the FT project
     * @param storeGoal       the goal providing the loaded database whose k-mer store is visited
     * @param bloomFilterGoal the goal providing the k-mer index Bloom filter used for
     *                        child-subtree membership tests
     * @param deps            further goals this goal depends on
     */
    @SafeVarargs
    public KMerBranchHistoGoal(P project, ObjectGoal<Database, P> storeGoal,
                               ObjectGoal<ProbFilter, P> bloomFilterGoal, Goal<P>... deps) {
        super(project, FTGoalKey.BRANCH_HISTO, storeGoal, bloomFilterGoal, deps);
    }

    /**
     * Initialises a fresh result map before the k-mer store traversal starts.
     */
    @Override
    protected void beforeKMerStoreWork() {
        map = new HashMap<>();
    }

    /**
     * Adds the current k-mer to its parent's branching-degree histogram.
     */
    @Override
    protected void inKMerStoreWork(SmallTaxIdNode parent, long pos, boolean[] bits, int spread,
                                   int[] setSlots) {
        int childCount = parent.getSubNodes().length;
        // The branching degree counts the parent's direct child subtrees plus the trailing "OTHER"
        // bucket (k-mers attributed to no child genome), which is treated like an additional child.
        // The base class' spread already includes the OTHER slot, so it is the branching degree
        // directly and can range from 1 up to childCount + 1.
        int degree = spread;
        if (degree >= 1) {
            long[] histo = map.get(parent.getTaxId());
            if (histo == null) {
                histo = new long[childCount + 1];
                map.put(parent.getTaxId(), histo);
            }
            histo[degree - 1]++;
        }
    }

    /**
     * Publishes the accumulated histogram map as this goal's result.
     */
    @Override
    protected void afterKMerStoreWork() {
        set(map);
        map = null;
    }
}
