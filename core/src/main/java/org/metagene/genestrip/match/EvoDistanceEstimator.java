package org.metagene.genestrip.match;

import it.unimi.dsi.fastutil.objects.Object2LongMap;
import org.metagene.genestrip.store.Database;
import org.metagene.genestrip.tax.SmallTaxTree;

import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

public abstract class EvoDistanceEstimator {
    public Map<SmallTaxTree.SmallTaxIdNode, DistanceInfo> computeDistances(Database database) {
        SmallTaxTree tree = database.getTaxTree();
        Object2LongMap<String> stats = database.getStats();

        Map<SmallTaxTree.SmallTaxIdNode, DistanceInfo> distances = new HashMap<>();
        Iterator<SmallTaxTree.SmallTaxIdNode> iterator = tree.iterator();
        while (iterator.hasNext()) {
            SmallTaxTree.SmallTaxIdNode node = iterator.next();
            distances.put(node, getDistanceForNode(node, stats));
        }

        return distances;
    }

    protected DistanceInfo getDistanceForNode(SmallTaxTree.SmallTaxIdNode node, Object2LongMap<String> stats) {
        DistanceInfo info = new DistanceInfo();
        long below = firstBelow(node, stats, info);
        long above = 0;
        node = node.getParent();
        while (node != null) {
            if (node.getTaxId() != null) {
                above += stats.getOrDefault(node.getTaxId(), 0L);
            }
            node = node.getParent();
        }
        long sum = above + below;
        info.distance = 1 - Math.pow(1 - ((double) below / sum), getInvK());
        return info;
    }

    protected abstract double getInvK();

    protected long firstBelow(SmallTaxTree.SmallTaxIdNode node, Object2LongMap<String> stats, DistanceInfo info) {
        long childMax = 0;
        SmallTaxTree.SmallTaxIdNode[] children = node.getSubNodes();
        if (children != null) {
            for (int i = 0; i < children.length; i++) {
                long down = below(children[i], stats);
                if (down > childMax) {
                    childMax = down;
                    info.branch = children[i];
                }
            }
        }
        return childMax + stats.getOrDefault(node.getTaxId(), 0L);
    }

    protected long below(SmallTaxTree.SmallTaxIdNode node, Object2LongMap<String> stats) {
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

    public static class DistanceInfo {
        private double distance;
        private SmallTaxTree.SmallTaxIdNode branch;

        public DistanceInfo() {
        }

        public double getDistance() {
            return distance;
        }

        public SmallTaxTree.SmallTaxIdNode getBranch() {
            return branch;
        }
    }
}
