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
package org.metagene.genestrip.finertree.cluster;

// This follows the simple standard algorithm "HAC" for agglomerative clustering
// as described in Figure 17.2. in Manning's "Introduction to Information Retrieval"
// Too inefficient for larger problems, but sufficient here.
/**
 * Performs hierarchical agglomerative clustering (HAC) over a set of items described by a
 * {@link Similarity} matrix. The implementation follows the simple standard algorithm from
 * Figure 17.2 of Manning's "Introduction to Information Retrieval": it repeatedly merges the two
 * most similar clusters until a single {@link DendrogramNode} tree remains. The way similarities
 * between a newly merged cluster and the remaining clusters are recomputed is controlled by the
 * configured {@link Method linkage method}. The algorithm runs in roughly quadratic time and space
 * and is therefore intended for small problems only.
 */
public class SimpleAggloClustering {
    /**
     * The linkage strategy used to recompute the similarity between a merged cluster and the
     * remaining clusters: nearest-neighbour ({@code SINGLE_LINKAGE}), farthest-neighbour
     * ({@code COMPLETE_LINKAGE}), size-weighted group average ({@code UPGMA}) or simple unweighted
     * group average ({@code WPGMA}).
     */
    public enum Method {
        /** Nearest-neighbour linkage: the similarity of the most similar pair of members. */
        SINGLE_LINKAGE,
        /** Farthest-neighbour linkage: the similarity of the least similar pair of members. */
        COMPLETE_LINKAGE,
        /** Group average weighted by cluster size. */
        UPGMA,
        /** Unweighted group average, i.e. both merged clusters count equally. */
        WPGMA
    };

    private final Method method;

    /**
     * Creates a clustering instance that merges clusters using the given linkage method.
     *
     * @param method the {@link Method linkage method} used to recompute cluster similarities
     */
    public SimpleAggloClustering(Method method) {
        this.method = method;
    }

    /**
     * Clusters all items described by the given similarity matrix and returns the root of the
     * resulting dendrogram. Each item initially forms its own leaf cluster; the two most similar
     * clusters are merged repeatedly until a single tree remains, with inter-cluster similarities
     * recomputed after every merge according to the configured {@link Method linkage method}.
     *
     * @param similarity the pairwise {@link Similarity} between the items to cluster
     * @return the root {@link DendrogramNode} of the fully merged dendrogram
     */
    public DendrogramNode cluster(Similarity similarity) {
        DendrogramNode[] clusters = new DendrogramNode[similarity.values()];
        // Only about half of this array is really needed. (Could optimize, but it's not worth it.)
        double[][] sims = new double[clusters.length][clusters.length];
        int[] sizes = new int[clusters.length];

        for (int i = 0; i < clusters.length; i++) {
            clusters[i] = new DendrogramNode(i, similarity.getSimilarity(i, i));
            sizes[i] = 1;
        }
        for (int i = 0; i < sims.length; i++) {
            for (int j = i + 1; j < sims.length; j++) {
                sims[j][i] = sims[i][j] = similarity.getSimilarity(i, j);
            }
        }

        double bestSim;
        int bestI = 0;
        int bestJ;
        for (int k = 0; k < clusters.length - 1; k++) {
            bestSim = -1; // Can't be zero, cause zero might occur as actual similarity value...
            bestI = 0;
            bestJ = 0;
            for (int i = 0; i < clusters.length; i++) {
                for (int j = i + 1; j < clusters.length; j++) {
                    if (clusters[i] != null && clusters[j] != null) {
                        if (sims[i][j] > bestSim) {
                            bestSim = sims[i][j];
                            bestI = i;
                            bestJ = j;
                        }
                    }
                }
            }
            DendrogramNode node = new DendrogramNode(clusters[bestI], clusters[bestJ], bestSim);
            clusters[bestI] = node;
            clusters[bestJ] = null;
            for (int h = 0; h < clusters.length; h++) {
                if (clusters[h] != null && h != bestI) {
                    sims[bestI][h] = sims[h][bestI] = similarity(similarity, sims, bestI, bestJ, h,sizes);
                }
            }
            sizes[bestI] += sizes[bestJ];
        }
        return clusters[bestI];
    }

    /**
     * Computes the updated similarity between the freshly merged cluster (the union of clusters
     * {@code bestI} and {@code bestJ}) and the remaining cluster {@code h}, dispatching to the
     * linkage helper that matches the configured {@link Method}.
     *
     * @param similarity the original pairwise {@link Similarity} of the items being clustered
     * @param sims       the current cluster-to-cluster similarity matrix
     * @param bestI      the index of the first cluster that was merged (and now holds the union)
     * @param bestJ      the index of the second cluster that was merged
     * @param h          the index of the remaining cluster whose similarity is being recomputed
     * @param sizes      the number of leaf items contained in each cluster
     * @return the updated similarity between the merged cluster and cluster {@code h}
     */
    protected double similarity(Similarity similarity, double[][] sims, int bestI, int bestJ, int h, int[] sizes) {
        switch (method) {
            case SINGLE_LINKAGE:
                return singleLinkage(sims, bestI, bestJ, h, sizes);
            case COMPLETE_LINKAGE:
                return completeLinkage(sims, bestI, bestJ, h, sizes);
            case UPGMA:
                return upgma(sims, bestI, bestJ, h, sizes);
            default:
                return wpgma(sims, bestI, bestJ, h, sizes);
        }
    }

    /**
     * Computes the single-linkage (nearest-neighbour) similarity between the merged cluster and
     * cluster {@code h} as the maximum of the two constituent clusters' similarities to {@code h}.
     *
     * @param sims  the current cluster-to-cluster similarity matrix
     * @param bestI the index of the first merged cluster
     * @param bestJ the index of the second merged cluster
     * @param h     the index of the remaining cluster
     * @param sizes the number of leaf items contained in each cluster
     * @return the nearest-neighbour similarity between the merged cluster and cluster {@code h}
     */
    protected double singleLinkage(double[][] sims, int bestI, int bestJ, int h, int[] sizes) {
        return Math.max(sims[h][bestI], sims[h][bestJ]);
    }

    /**
     * Computes the complete-linkage (farthest-neighbour) similarity between the merged cluster and
     * cluster {@code h} as the minimum of the two constituent clusters' similarities to {@code h}.
     *
     * @param sims  the current cluster-to-cluster similarity matrix
     * @param bestI the index of the first merged cluster
     * @param bestJ the index of the second merged cluster
     * @param h     the index of the remaining cluster
     * @param sizes the number of leaf items contained in each cluster
     * @return the farthest-neighbour similarity between the merged cluster and cluster {@code h}
     */
    protected double completeLinkage(double[][] sims, int bestI, int bestJ, int h, int[] sizes) {
        return Math.min(sims[h][bestI], sims[h][bestJ]);
    }

    // According to:
    // https://en.wikipedia.org/wiki/UPGMA
    // Seems to the same as "group average"
    /**
     * Computes the UPGMA (unweighted pair group method with arithmetic mean, i.e. group-average)
     * similarity between the merged cluster and cluster {@code h} as the size-weighted average of
     * the two constituent clusters' similarities to {@code h}.
     *
     * @param sims  the current cluster-to-cluster similarity matrix
     * @param bestI the index of the first merged cluster
     * @param bestJ the index of the second merged cluster
     * @param h     the index of the remaining cluster
     * @param sizes the number of leaf items contained in each cluster, used as averaging weights
     * @return the group-average similarity between the merged cluster and cluster {@code h}
     */
    protected double upgma(double[][] sims, int bestI, int bestJ, int h, int[] sizes) {
        return (sizes[bestI] * sims[h][bestI] + sizes[bestJ] * sims[h][bestJ]) / (sizes[bestI] + sizes[bestJ]);
    }

    /**
     * Computes the WPGMA (weighted pair group method with arithmetic mean) similarity between the
     * merged cluster and cluster {@code h} as the simple, size-independent mean of the two
     * constituent clusters' similarities to {@code h}.
     *
     * @param sims  the current cluster-to-cluster similarity matrix
     * @param bestI the index of the first merged cluster
     * @param bestJ the index of the second merged cluster
     * @param h     the index of the remaining cluster
     * @param sizes the number of leaf items contained in each cluster (unused by this method)
     * @return the unweighted average similarity between the merged cluster and cluster {@code h}
     */
    protected double wpgma(double[][] sims, int bestI, int bestJ, int h, int[] sizes) {
        return (sims[h][bestI] + sims[h][bestJ]) / 2;
    }
}
