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
package org.metagene.genestrip.match;

import it.unimi.dsi.fastutil.objects.Object2LongMap;
import org.metagene.genestrip.store.Database;
import org.metagene.genestrip.tax.SmallTaxTree;

import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

public class EvoDistanceEstimator {
    public Map<SmallTaxTree.SmallTaxIdNode, DistanceInfo> computeDistances(Database database) {
        SmallTaxTree tree = database.getTaxTree();
        Object2LongMap<String> stats = database.getStats();

        Map<SmallTaxTree.SmallTaxIdNode, DistanceInfo> distances = new HashMap<>();
        Iterator<SmallTaxTree.SmallTaxIdNode> iterator = tree.iterator();
        int k = database.getKmerStore().getK();
        while (iterator.hasNext()) {
            SmallTaxTree.SmallTaxIdNode node = iterator.next();
            distances.put(node, new DistanceInfo(node, stats, k));
        }
        for (DistanceInfo distanceInfo : distances.values()) {
            distanceInfo.updateDistancePortion(distances.get(distanceInfo.getBranch()));
        }

        return distances;
    }

    public static class DistanceInfo {
        private final SmallTaxTree.SmallTaxIdNode node;
        private final double distance;
        private final SmallTaxTree.SmallTaxIdNode branch;
        private double distancePortion;

        public DistanceInfo(SmallTaxTree.SmallTaxIdNode node, Object2LongMap<String> stats, int k) {
            this.node = node;

            long childMax = 0;
            SmallTaxTree.SmallTaxIdNode bestBranch = null;
            SmallTaxTree.SmallTaxIdNode[] children = node.getSubNodes();
            if (children != null) {
                for (int i = 0; i < children.length; i++) {
                    long down = below(children[i], stats);
                    if (down > childMax) {
                        childMax = down;
                        bestBranch = children[i];
                    }
                }
            }
            branch = bestBranch;

            long below = childMax + stats.getOrDefault(node.getTaxId(), 0L);
            long above = 0;
            node = node.getParent();
            while (node != null) {
                if (node.getTaxId() != null) {
                    above += stats.getOrDefault(node.getTaxId(), 0L);
                }
                node = node.getParent();
            }
            long sum = above + below;
            distance = 1 - Math.pow(1 - ((double) below / sum), 1d / k);
        }

        public SmallTaxTree.SmallTaxIdNode getNode() {
            return node;
        }

        private void updateDistancePortion(DistanceInfo childInfo) {
            distancePortion = distance - (childInfo == null ? 0 : childInfo.distance);
        }

        private long below(SmallTaxTree.SmallTaxIdNode node, Object2LongMap<String> stats) {
            long childMax = 0;
            if (node.getSubNodes() != null) {
                for (SmallTaxTree.SmallTaxIdNode subNode : node.getSubNodes()) {
                    long down = below(subNode, stats);
                    if (down > childMax) {
                        childMax = down;
                    }
                }
            }
            return childMax + stats.getOrDefault(node.getTaxId(), 0L);
        }

        public double getDistance() {
            return distance;
        }

        public double getDistancePortion() {
            return distancePortion;
        }

        public SmallTaxTree.SmallTaxIdNode getBranch() {
            return branch;
        }
    }
}
