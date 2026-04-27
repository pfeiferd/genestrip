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
public class SimpleAggloClustering {
    public enum Method { SINGLE_LINKAGE, COMPLETE_LINKAGE, UPGMA, WPGMA};

    private final Method method;

    public SimpleAggloClustering(Method method) {
        this.method = method;
    }

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
            sizes[bestI] += sizes[bestJ];
            for (int h = 0; h < clusters.length; h++) {
                if (clusters[h] != null && h != bestI) {
                    sims[bestI][h] = sims[h][bestI] = similarity(similarity, sims, bestI, bestJ, h,sizes);
                }
            }
        }
        return clusters[bestI];
    }

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

    protected double singleLinkage(double[][] sims, int bestI, int bestJ, int h, int[] sizes) {
        return Math.max(sims[h][bestI], sims[h][bestJ]);
    }

    protected double completeLinkage(double[][] sims, int bestI, int bestJ, int h, int[] sizes) {
        return Math.min(sims[h][bestI], sims[h][bestJ]);
    }

    // According to:
    // https://en.wikipedia.org/wiki/UPGMA
    // Seems to the same as "group average"
    protected double upgma(double[][] sims, int bestI, int bestJ, int h, int[] sizes) {
        return (sizes[bestI] * sims[h][bestI] + sizes[bestJ] * sims[h][bestJ]) / (sizes[bestI] + sizes[bestJ]);
    }

    protected double wpgma(double[][] sims, int bestI, int bestJ, int h, int[] sizes) {
        return (sims[h][bestI] + sims[h][bestJ]) / 2;
    }
}
