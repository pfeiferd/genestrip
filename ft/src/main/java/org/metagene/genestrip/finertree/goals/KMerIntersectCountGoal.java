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
import org.metagene.genestrip.finertree.FTConfigKey;
import org.metagene.genestrip.finertree.FTGoalKey;
import org.metagene.genestrip.finertree.FTProject;
import org.metagene.genestrip.finertree.bloom.XORKMerIndexBloomFilter;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.store.Database;
import org.metagene.genestrip.tax.SmallTaxTree;

import java.util.*;

public class KMerIntersectCountGoal<P extends FTProject> extends KMerStoreWorkGoal<KMerIntersectCountGoal.IntersectionsPerNode, P> {
    public interface IntersectionsPerNode  {
        public Set<SmallTaxTree.SmallTaxIdNode> getParentNodes();
        public long getIntersectionCount(SmallTaxTree.SmallTaxIdNode parent, int child1, int child2);
        public long getKMerSpreadSum(SmallTaxTree.SmallTaxIdNode parent);
        public long getKMerSum(SmallTaxTree.SmallTaxIdNode parent);
        public double getJaccardIndex(SmallTaxTree.SmallTaxIdNode parent, int i, int j, boolean withChildCounts);
        public double getAvgKMerSpread(SmallTaxTree.SmallTaxIdNode parent);
        public double getOverspreadRatio(SmallTaxTree.SmallTaxIdNode parent);
        public long getSubnodesKMerCount(SmallTaxTree.SmallTaxIdNode parent);
    }

    private IntersectionsPerNodeImpl intersectionsPerNode;

    @SafeVarargs
    public KMerIntersectCountGoal(P project, ObjectGoal<Database, P> storeGoal,
                              ObjectGoal<XORKMerIndexBloomFilter, P> bloomFilterGoal,
                              Goal<P>... deps) {
        super(project, FTGoalKey.INTERSECT_COUNT, storeGoal, bloomFilterGoal, deps);
    }

    @Override
    protected void beforeKMerStoreWork() {
        intersectionsPerNode = new IntersectionsPerNodeImpl();
    }

    @Override
    protected void inKMerStoreWork(SmallTaxTree.SmallTaxIdNode parent, long pos, boolean[] bits, int spread) {
        SmallTaxTree.SmallTaxIdNode[] children = parent.getSubNodes();
        intersectionsPerNode.incKMerSpread(parent, spread);
        for (int i = 0; i <= children.length; i++) {
            for (int j = i; j <= children.length; j++) {
                if (bits[i] && bits[j]) {
                    intersectionsPerNode.incIntersectionCount(parent, i, j);
                }
            }
        }
    }

    @Override
    protected void afterKMerStoreWork() {
        Object2LongMap<SmallTaxTree.SmallTaxIdNode> stats = kMerSortedArray.getNKmersPerTaxid();
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

    public class IntersectionsPerNodeImpl implements IntersectionsPerNode {
        private Set<SmallTaxTree.SmallTaxIdNode> immutableParentNodes;
        private Map<SmallTaxTree.SmallTaxIdNode, long[]> parentToCounts;
        private Map<SmallTaxTree.SmallTaxIdNode, long[]> childToSubnodeCounts;

        public IntersectionsPerNodeImpl() {
            parentToCounts = new HashMap<>();
            childToSubnodeCounts = new HashMap<>();
            immutableParentNodes = Collections.unmodifiableSet(parentToCounts.keySet());
        }

        @Override
        public Set<SmallTaxTree.SmallTaxIdNode> getParentNodes() {
            return immutableParentNodes;
        }

        @Override
        public long getIntersectionCount(SmallTaxTree.SmallTaxIdNode parent, int i, int j) {
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

        private long[] countsForParent(SmallTaxTree.SmallTaxIdNode parent) {
            long[] counts = parentToCounts.get(parent);
            if (counts == null) {
                // "+ 1" for "OTHER_VALUE", spread and number of k-mers
                int c = parent.getSubNodes().length + 1;
                parentToCounts.put(parent, counts = new long[(c * c + c) / 2 + 2]);
            }
            return counts;
        }

        void incKMerSpread(SmallTaxTree.SmallTaxIdNode parent, int spread) {
            long[] counts = countsForParent(parent);
            counts[counts.length - 2] += spread;
            counts[counts.length - 1]++;
        }

        // No Laplace correction
        @Override
        public long getKMerSpreadSum(SmallTaxTree.SmallTaxIdNode parent) {
            long[] counts = parentToCounts.get(parent);
            return counts == null ? 0 : counts[counts.length - 2];
        }

        // No Laplace correction
        @Override
        public long getKMerSum(SmallTaxTree.SmallTaxIdNode parent) {
            long[] counts = parentToCounts.get(parent);
            return counts == null ? 0 : counts[counts.length - 1];
        }

        // With Laplace correction
        @Override
        public double getJaccardIndex(SmallTaxTree.SmallTaxIdNode parent, int i, int j, boolean withChildCounts) {
            long intersect = getIntersectionCount(parent, i, j) + 1; // "+ 1" is Laplace smoothing
            long union = getIntersectionCount(parent, i, i) + 1 + getIntersectionCount(parent, j, j) + 1 - intersect; // "+ 1" is Laplace smoothing
            if (withChildCounts) {
                SmallTaxTree.SmallTaxIdNode[] children = parent.getSubNodes();
                if (i < children.length && j < children.length) {
                    union += getSubnodesKMerCount(children[i]);
                    union += getSubnodesKMerCount(children[j]);
                }
            }
            if (union == 0) {
                if (intersect != 0) {
                    throw new IllegalStateException("If union is 0 then intersect must be too.");
                }
                // Both sets are empty - so they are identical and perfectly similar...
                return 1;
            }
            return ((double) intersect) / union;
        }

        @Override
        public double getAvgKMerSpread(SmallTaxTree.SmallTaxIdNode parent) {
            return ((double) getKMerSpreadSum(parent)) / getKMerSum(parent);
        }

        @Override
        public double getOverspreadRatio(SmallTaxTree.SmallTaxIdNode parent) {
            return (getAvgKMerSpread(parent) - 2) / ((parent.getSubNodes().length + 1) - 2);
        }

        @Override
        public long getSubnodesKMerCount(SmallTaxTree.SmallTaxIdNode child) {
            long[] count = childToSubnodeCounts.get(child);
            return count == null ? 0 : count[0];
        }
    }
}
