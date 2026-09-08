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
import org.metagene.genestrip.probfilter.ProbFilter;
import org.metagene.genestrip.finertree.FTConfigKey;
import org.metagene.genestrip.finertree.FTGoalKey;
import org.metagene.genestrip.finertree.FTProject;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.store.Database;
import org.metagene.genestrip.tax.SmallTaxTree;

import java.util.Collections;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

/**
 * Concrete {@link KMerStoreWorkGoal} that, per taxonomy parent node at a refinement rank, accumulates
 * pairwise intersection counts of *k*-mer presence across the parent's child subtrees together with
 * *k*-mer spread statistics. For every visited *k*-mer it increments the intersection count for each
 * pair of set child-membership slots (including a child paired with itself and with the trailing
 * "OTHER" slot) and records the *k*-mer's spread. The accumulated result is exposed as an
 * {@link IntersectionsPerNode} instance from which intersection, union, Jaccard and spread values can
 * be derived.
 *
 * @param <P> the concrete FT project type
 */
public class KMerIntersectCountGoal<P extends FTProject> extends KMerStoreWorkGoal<KMerIntersectCountGoal.IntersectionsPerNode, P> {
    /**
     * Read-only view of the intersection and spread statistics accumulated per parent node. Child
     * subtrees are addressed by index; the index equal to the number of children denotes the trailing
     * "OTHER" slot.
     */
    public interface IntersectionsPerNode  {
        /**
         * Returns the set of parent nodes for which statistics have been collected.
         *
         * @return the parent nodes
         */
        public Set<SmallTaxTree.SmallTaxIdNode> getParentNodes();

        /**
         * Returns the number of *k*-mers present in both of the given child subtrees of the parent
         * (or, when {@code child1 == child2}, the number of *k*-mers present in that single subtree).
         *
         * <p>
         * When {@code withDescendantCounts} is set and {@code child1 == child2}, the descendant *k*-mer
         * counts of that child subtree (recursively accumulated in {@link #getSubnodesKMerCount}) are
         * added to the count. Off-diagonal pairs and the "OTHER" slot are unaffected.
         *
         * @param parent               the parent node
         * @param child1               the index of the first child subtree (or "OTHER" slot)
         * @param child2               the index of the second child subtree (or "OTHER" slot)
         * @param withDescendantCounts whether to add the child's descendant *k*-mer counts on the diagonal
         * @return the intersection count, or {@code 0} if the parent is unknown
         */
        public long getIntersectionCount(SmallTaxTree.SmallTaxIdNode parent, int child1, int child2, boolean withDescendantCounts);

        /**
         * Returns the sum of the spread values (number of subtrees each *k*-mer occurs in) over all
         * *k*-mers of the parent, without Laplace correction.
         *
         * @param parent the parent node
         * @return the total *k*-mer spread, or {@code 0} if the parent is unknown
         */
        public long getKMerSpreadSum(SmallTaxTree.SmallTaxIdNode parent);

        /**
         * Returns the number of *k*-mers counted for the parent, without Laplace correction.
         *
         * @param parent the parent node
         * @return the *k*-mer count, or {@code 0} if the parent is unknown
         */
        public long getKMerSum(SmallTaxTree.SmallTaxIdNode parent);

        /**
         * Returns the Laplace-smoothed Jaccard index between the two given child subtrees of the parent,
         * optionally enlarging the two set sizes by the descendant *k*-mer counts of the two children.
         *
         * @param parent               the parent node
         * @param i                    the index of the first child subtree (or "OTHER" slot)
         * @param j                    the index of the second child subtree (or "OTHER" slot)
         * @param withDescendantCounts whether to add the children's descendant *k*-mer counts to the union
         * @return the Jaccard index in {@code [0, 1]}; {@code 1} when both sets are empty
         */
        public double getJaccardIndex(SmallTaxTree.SmallTaxIdNode parent, int i, int j, boolean withDescendantCounts);

        /**
         * Returns the average *k*-mer spread for the parent, i.e. the spread sum divided by the *k*-mer
         * count.
         *
         * @param parent the parent node
         * @return the average *k*-mer spread
         */
        public double getAvgKMerSpread(SmallTaxTree.SmallTaxIdNode parent);

        /**
         * Returns the overspread ratio for the parent, normalising the average spread above the minimum
         * of two to the range spanned up to the maximum possible spread (number of children plus the
         * "OTHER" slot).
         *
         * @param parent the parent node
         * @return the overspread ratio
         */
        public double getOverspreadRatio(SmallTaxTree.SmallTaxIdNode parent);

        /**
         * Returns the number of *k*-mers assigned to the entire subtree rooted at the given node.
         *
         * @param parent the node whose subtree *k*-mer count is requested
         * @return the subtree *k*-mer count, or {@code 0} if unknown
         */
        public long getSubnodesKMerCount(SmallTaxTree.SmallTaxIdNode parent);
    }

    private IntersectionsPerNodeImpl intersectionsPerNode;
    /** Whether to merge on the Jaccard index rather than containment, see
     * {@link FTConfigKey#JACCARD_SIM}. */
    private final boolean jaccardSim;

    /**
     * Creates the goal under the {@link FTGoalKey#INTERSECT_COUNT} key.
     *
     * @param project         the FT project
     * @param storeGoal       the goal providing the loaded database whose *k*-mer store is visited
     * @param bloomFilterGoal the goal providing the XOR *k*-mer index bloom filter
     * @param deps            further goals this goal depends on
     */
    @SafeVarargs
    public KMerIntersectCountGoal(P project, ObjectGoal<Database, P> storeGoal,
                              ObjectGoal<ProbFilter, P> bloomFilterGoal,
                              Goal<P>... deps) {
        super(project, FTGoalKey.INTERSECT_COUNT, storeGoal, bloomFilterGoal, deps);
        jaccardSim = booleanConfigValue(FTConfigKey.JACCARD_SIM);
    }

    /**
     * Initialises a fresh {@link IntersectionsPerNodeImpl} accumulator before the traversal starts.
     */
    @Override
    protected void beforeKMerStoreWork() {
        intersectionsPerNode = new IntersectionsPerNodeImpl();
    }

    /**
     * Records the *k*-mer's spread for the parent and increments the intersection count for every pair
     * of set child-membership slots (including each slot paired with itself and with the "OTHER" slot).
     */
    @Override
    protected void inKMerStoreWork(SmallTaxTree.SmallTaxIdNode parent, long pos, boolean[] bits, int spread,
                                   int[] setSlots) {
        // Looked up once for this k-mer and then indexed directly. It used to be fetched from a map
        // inside the loop below - once per counted pair, plus once for the spread.
        long[] counts = intersectionsPerNode.countsForParent(parent);
        counts[counts.length - 2] += spread;
        counts[counts.length - 1]++;
        // Over the slots that are actually set, and not over every pair of slots the node has. The
        // pairs a k-mer contributes to are those among the children carrying it, and a k-mer is
        // typically carried by a handful of them: the old loop ran (c+1)^2/2 times whatever the
        // spread, which for a node with 3,500 children is six million iterations per k-mer to count
        // a few. Ascending slots, so i <= j holds and the triangular index needs no ordering.
        for (int a = 0; a < spread; a++) {
            int i = setSlots[a];
            for (int b = a; b < spread; b++) {
                int j = setSlots[b];
                counts[(j * j + j) / 2 + i]++;
            }
        }
    }

    /**
     * Distributes the store's per-tax-id *k*-mer counts to the nearest ancestor parent that sits at a
     * refinement rank, thereby populating the descendant *k*-mer counts, and finally publishes the
     * accumulated {@link IntersectionsPerNode} as this goal's result.
     */
    @Override
    protected void afterKMerStoreWork() {
        Object2LongMap<SmallTaxTree.SmallTaxIdNode> stats = kMerStore.getNKmersPerTaxid();
        stats.forEach((s, aLong) -> {
            while (s != null) {
                SmallTaxTree.SmallTaxIdNode parent = s.getParent();
                if (parent != null) {
                    if (FTConfigKey.RefinementPosition.getMatchingNodeFor(parent, refinementIntervals) != null) {
                        intersectionsPerNode.incSubnodeCounts(s, aLong);
                        break;
                    }
                }
                s = parent;
            }
        });

        set(intersectionsPerNode);
    }

    /**
     * Default {@link IntersectionsPerNode} implementation backed by per-parent {@code long[]} arrays.
     * Pairwise intersection counts are stored in a triangular layout indexed by {@code (j*j+j)/2 + i}
     * (with {@code i <= j}); the two trailing array entries hold the accumulated *k*-mer spread sum and
     * *k*-mer count. A separate map records the descendant *k*-mer count for each child subtree.
     */
    public class IntersectionsPerNodeImpl implements IntersectionsPerNode {
        private Set<SmallTaxTree.SmallTaxIdNode> immutableParentNodes;
        private Map<SmallTaxTree.SmallTaxIdNode, long[]> parentToCounts;
        private Map<SmallTaxTree.SmallTaxIdNode, long[]> childToSubnodeCounts;

        /**
         * Creates an empty accumulator with no parent statistics collected yet.
         */
        public IntersectionsPerNodeImpl() {
            parentToCounts = new HashMap<>();
            childToSubnodeCounts = new HashMap<>();
            immutableParentNodes = Collections.unmodifiableSet(parentToCounts.keySet());
        }

        /**
         * Returns an unmodifiable view of the parent nodes for which counts have been accumulated.
         */
        @Override
        public Set<SmallTaxTree.SmallTaxIdNode> getParentNodes() {
            return immutableParentNodes;
        }

        /**
         * Looks up the intersection count for the pair of slots {@code (i, j)} and, when
         * {@code withDescendantCounts} is set on the diagonal ({@code i == j}), adds the child's
         * recursively accumulated descendant *k*-mer counts. There is no descendant data for the cross
         * intersection of two distinct children, so off-diagonal pairs and the "OTHER" slot are returned
         * unchanged.
         */
        @Override
        public long getIntersectionCount(SmallTaxTree.SmallTaxIdNode parent, int i, int j, boolean withDescendantCounts) {
            long count = rawIntersectionCount(parent, i, j);
            if (withDescendantCounts && i == j) {
                SmallTaxTree.SmallTaxIdNode[] children = parent.getSubNodesWithoutOther();
                if (i < children.length) {
                    count += getSubnodesKMerCount(children[i]);
                }
            }
            return count;
        }

        /**
         * Looks up the raw intersection count for the pair of slots {@code (i, j)} in the parent's
         * triangular count array, normalising the order of the indices first and returning {@code 0}
         * when no counts exist for the parent.
         */
        private long rawIntersectionCount(SmallTaxTree.SmallTaxIdNode parent, int i, int j) {
            if (i > j) {
                int h = i;
                i = j;
                j = h;
            }
            long[] counts = parentToCounts.get(parent);
            return counts == null ? 0 : counts[(j * j + j) / 2 + i];
        }

        void incSubnodeCounts(SmallTaxTree.SmallTaxIdNode parent, long add) {
            long[] count = childToSubnodeCounts.get(parent);
            if (count == null) {
                count = new long[1];
                childToSubnodeCounts.put(parent, count);
            }
            count[0] += add;
        }

        void incIntersectionCount(SmallTaxTree.SmallTaxIdNode parent, int i, int j) {
            long[] counts = countsForParent(parent);
            if (i > j) {
                int h = i;
                i = j;
                j = h;
            }
            counts[(j * j + j) / 2 + i]++;
        }

        long[] countsForParent(SmallTaxTree.SmallTaxIdNode parent) {
            long[] counts = parentToCounts.get(parent);
            if (counts == null) {
                // "+ 1" for "OTHER_VALUE", spread and number of k-mers
                int c = parent.getSubNodesWithoutOther().length + 1;
                parentToCounts.put(parent, counts = new long[(c * c + c) / 2 + 2]);
            }
            return counts;
        }

        void incKMerSpread(SmallTaxTree.SmallTaxIdNode parent, int spread) {
            long[] counts = countsForParent(parent);
            counts[counts.length - 2] += spread;
            counts[counts.length - 1]++;
        }

        /**
         * Returns the accumulated *k*-mer spread sum stored in the second-to-last entry of the parent's
         * count array, or {@code 0} when the parent is unknown. No Laplace correction is applied.
         */
        // No Laplace correction
        @Override
        public long getKMerSpreadSum(SmallTaxTree.SmallTaxIdNode parent) {
            long[] counts = parentToCounts.get(parent);
            return counts == null ? 0 : counts[counts.length - 2];
        }

        /**
         * Returns the accumulated *k*-mer count stored in the last entry of the parent's count array,
         * or {@code 0} when the parent is unknown. No Laplace correction is applied.
         */
        // No Laplace correction
        @Override
        public long getKMerSum(SmallTaxTree.SmallTaxIdNode parent) {
            long[] counts = parentToCounts.get(parent);
            return counts == null ? 0 : counts[counts.length - 1];
        }

        /**
         * Computes the Jaccard index of slots {@code i} and {@code j} as {@code |A ∩ B| / |A ∪ B|} using
         * the stored intersection counts, adding one to intersection and to each self-count for Laplace
         * smoothing. When {@code withDescendantCounts} is set, the recursively accumulated descendant
         * *k*-mer counts extend the two set sizes via the diagonal
         * {@link #getIntersectionCount(SmallTaxTree.SmallTaxIdNode, int, int, boolean) intersection counts}.
         * Returns {@code 1} when both sets are empty.
         */
        // With Laplace correction
        @Override
        public double getJaccardIndex(SmallTaxTree.SmallTaxIdNode parent, int i, int j, boolean withDescendantCounts) {
            long intersect = getIntersectionCount(parent, i, j, withDescendantCounts) + 1; // "+ 1" is Laplace smoothing
            long sizeI = getIntersectionCount(parent, i, i, withDescendantCounts) + 1; // "+ 1" is Laplace smoothing
            long sizeJ = getIntersectionCount(parent, j, j, withDescendantCounts) + 1; // "+ 1" is Laplace smoothing
            // Both denominators are >= 1 after Laplace smoothing, so no zero-division guard is needed.
            if (!jaccardSim) {
                // How much of the smaller of the two sets the two share. A child whose k-mer set is
                // small because its assembly is fragmented - sequence missing rather than differing -
                // is then not pushed away from everything it belongs with, which dividing by the union
                // does: its union with any well-assembled genome is nearly that genome's whole set.
                return ((double) intersect) / Math.min(sizeI, sizeJ);
            }
            return ((double) intersect) / (sizeI + sizeJ - intersect);
        }

        /**
         * Returns the spread sum divided by the *k*-mer count for the parent.
         */
        @Override
        public double getAvgKMerSpread(SmallTaxTree.SmallTaxIdNode parent) {
            return ((double) getKMerSpreadSum(parent)) / getKMerSum(parent);
        }

        /**
         * Returns the average spread rescaled so that the minimum spread of two maps to {@code 0} and the
         * maximum possible spread (number of children plus the "OTHER" slot) maps to {@code 1}.
         */
        @Override
        public double getOverspreadRatio(SmallTaxTree.SmallTaxIdNode parent) {
            return (getAvgKMerSpread(parent) - 2) / ((parent.getSubNodesWithoutOther().length + 1) - 2);
        }

        /**
         * Returns the descendant *k*-mer count recorded for the given child node, or {@code 0} when none
         * was recorded.
         */
        @Override
        public long getSubnodesKMerCount(SmallTaxTree.SmallTaxIdNode child) {
            long[] count = childToSubnodeCounts.get(child);
            return count == null ? 0 : count[0];
        }
    }
}
