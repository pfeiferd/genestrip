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

import org.junit.Test;
import org.metagene.genestrip.tax.Rank;
import org.metagene.genestrip.tax.SmallTaxTree;
import org.metagene.genestrip.tax.TaxTree;

import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.Map;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

/**
 * Tests the aggregation behind the subtree-relative quality measures, i.e.
 * {@link DBQualityCountsGoal#aggregateCounts} together with
 * {@link DBQualityCountsGoal.Counts#getSubtreePrecision()} and
 * {@link DBQualityCountsGoal.Counts#getPrecision()}.
 *
 * The base tree is
 * <pre>
 *   1 (root)
 *   └── 2
 *       ├── 3
 *       │   ├── 5 (DATA)
 *       │   └── 6 (DATA)
 *       └── 4 (DATA)
 * </pre>
 * with the k-mer assignment set up in {@link #baseCounts()}. The refined variant additionally inserts
 * the node 7 between 3 and its children - an added node, which the measure must not reward as long as
 * no k-mer actually moves.
 */
public class DBQualityCountsGoalTest {
    /** Tax ids of the data taxa, i.e. of the nodes of rank {@link Rank#DATA}. */
    private static final String[] DATA_TAX_IDS = {"4", "5", "6"};

    /**
     * Builds a {@link SmallTaxTree} from the given child-to-parent edges, giving the nodes listed in
     * {@link #DATA_TAX_IDS} the rank {@link Rank#DATA} and all others no rank.
     *
     * @param edges the tree edges as {@code {child, parent}} pairs, the root being its own parent
     * @return the resulting small tax tree with all nodes marked as required
     */
    private SmallTaxTree buildTree(int[][] edges) throws IOException {
        File dir = Files.createTempDirectory("dbqualitycounts").toFile();
        dir.deleteOnExit();
        StringBuilder nodes = new StringBuilder();
        StringBuilder names = new StringBuilder();
        for (int[] e : edges) {
            String rank = isDataTaxon(String.valueOf(e[0])) ? Rank.DATA.getName() : "no rank";
            nodes.append(e[0]).append("\t|\t").append(e[1]).append("\t|\t").append(rank).append("\t|\t\t|\n");
            names.append(e[0]).append("\t|\t").append(e[0]).append("\t|\t\t|\tscientific name\t|\n");
        }
        Files.write(new File(dir, TaxTree.NODES_DMP).toPath(), nodes.toString().getBytes(StandardCharsets.UTF_8));
        Files.write(new File(dir, TaxTree.NAMES_DMP).toPath(), names.toString().getBytes(StandardCharsets.UTF_8));
        TaxTree full = new TaxTree(dir, false);
        for (int[] e : edges) {
            full.getNodeByTaxId(String.valueOf(e[0])).markRequired();
        }
        return full.toSmallTaxTree();
    }

    private static boolean isDataTaxon(String taxId) {
        for (String dataTaxId : DATA_TAX_IDS) {
            if (dataTaxId.equals(taxId)) {
                return true;
            }
        }
        return false;
    }

    /**
     * Creates a tally for a node, seeded with the number of k-mers stored at it and with the number of
     * (k-mer, data taxon) pairs where the k-mer occurs in the taxon's genome.
     *
     * @param taxId    the tax id of the node
     * @param kmers    the number of k-mers stored at the node, i.e. |A_n|
     * @param tpOfNode the number of supporting (k-mer, data taxon) pairs, i.e. the sum of c_n(a)
     * @return the seeded tally
     */
    private DBQualityCountsGoal.Counts counts(String taxId, long kmers, long tpOfNode) {
        DBQualityCountsGoal.Counts counts = new DBQualityCountsGoal.Counts(isDataTaxon(taxId), kmers);
        for (long i = 0; i < tpOfNode; i++) {
            counts.incTpForNodePrecision();
        }
        return counts;
    }

    /**
     * The k-mer assignment shared by both tree variants. The root holds four k-mers occurring in all
     * three data taxa; they lie above node 2 and must not influence its weighted path precision. Node 6
     * holds no k-mers at all.
     */
    private Map<String, DBQualityCountsGoal.Counts> baseCounts() {
        Map<String, DBQualityCountsGoal.Counts> map = new LinkedHashMap<>();
        map.put("1", counts("1", 4, 12));
        map.put("2", counts("2", 2, 4));
        map.put("3", counts("3", 3, 6));
        map.put("4", counts("4", 2, 2));
        map.put("5", counts("5", 1, 1));
        map.put("6", counts("6", 0, 0));
        return map;
    }

    @Test
    public void testWeightedPathPrecisionIsRelativeToItsNode() throws IOException {
        SmallTaxTree tree = buildTree(new int[][]{{1, 1}, {2, 1}, {3, 2}, {4, 2}, {5, 3}, {6, 3}});
        Map<String, DBQualityCountsGoal.Counts> map = baseCounts();
        DBQualityCountsGoal.aggregateCounts(tree, map);

        // Data taxa underneath: |D_2| = 3, |D_3| = 2, |D_4| = |D_5| = |D_6| = 1.
        DBQualityCountsGoal.Counts m = map.get("2");
        // Pooled over the subtree of 2 only: 3*2 + 2*3 + 1*2 + 1*1 + 1*0 = 15, none of it from the root.
        assertEquals(15, m.getSubtreeTpPlusFp());
        assertEquals(4 + 6 + 2 + 1 + 0, m.getSubtreeTp());
        assertEquals(13.0 / 15.0, m.getPrecision(), 1e-12);

        // The root additionally pools its own four k-mers, shared by all three data taxa.
        DBQualityCountsGoal.Counts root = map.get("1");
        assertEquals(15 + 3 * 4, root.getSubtreeTpPlusFp());
        assertEquals(13 + 12, root.getSubtreeTp());
        assertEquals(25.0 / 27.0, root.getPrecision(), 1e-12);

        // The data taxa counted underneath each node, i.e. |D_n|.
        assertEquals(3, root.getLeaves());
        assertEquals(3, m.getLeaves());
        assertEquals(2, map.get("3").getLeaves());
        assertEquals(1, map.get("6").getLeaves());
    }

    @Test
    public void testSubtreePrecisionAveragesOverKMers() throws IOException {
        SmallTaxTree tree = buildTree(new int[][]{{1, 1}, {2, 1}, {3, 2}, {4, 2}, {5, 3}, {6, 3}});
        Map<String, DBQualityCountsGoal.Counts> map = baseCounts();
        DBQualityCountsGoal.aggregateCounts(tree, map);

        // Two k-mers at node 2, each occurring in two of its three data taxa.
        assertEquals(4.0 / 6.0, map.get("2").getNodePrecision(), 1e-12);
        // Node 6 holds no k-mers at all, so its node precision is undefined rather than perfect.
        assertTrue(Double.isNaN(map.get("6").getNodePrecision()));

        // Equivalently to averaging p(a) over the subtree's k-mers, the subtree precision is the mean
        // of the node precisions weighted by each node's share of those k-mers; nodes without k-mers
        // take part in neither the mean nor the shares.
        double weightedSum = 0;
        long kmerSum = 0;
        for (String taxId : new String[]{"2", "3", "4", "5", "6"}) {
            DBQualityCountsGoal.Counts counts = map.get(taxId);
            if (counts.getKmerSumForNode() > 0) {
                weightedSum += counts.getKmerSumForNode() * counts.getNodePrecision();
                kmerSum += counts.getKmerSumForNode();
            }
        }
        assertEquals(2 + 3 + 2 + 1, kmerSum);
        assertEquals(kmerSum, map.get("2").getSubtreeKmerSum());
        assertEquals(weightedSum / kmerSum, map.get("2").getSubtreePrecision(), 1e-12);

        // The definition proper: the plain mean of p(a) = c(a) / |D| over the subtree's k-mers, where
        // node 2 contributes 2 k-mers at 2/3, node 3 three at 1, node 4 two at 1 and node 5 one at 1.
        assertEquals((2 * (2.0 / 3.0) + 3 + 2 + 1) / 8.0, map.get("2").getSubtreePrecision(), 1e-12);

        // The root's own four k-mers enter its own average but never that of node 2.
        assertEquals(8 + 4, map.get("1").getSubtreeKmerSum());
        assertEquals((2 * (2.0 / 3.0) + 3 + 2 + 1 + 4) / 12.0, map.get("1").getSubtreePrecision(), 1e-12);
    }

    @Test
    public void testAddedNodesDoNotRaiseWeightedPathPrecision() throws IOException {
        SmallTaxTree plain = buildTree(new int[][]{{1, 1}, {2, 1}, {3, 2}, {4, 2}, {5, 3}, {6, 3}});
        Map<String, DBQualityCountsGoal.Counts> plainMap = baseCounts();
        DBQualityCountsGoal.aggregateCounts(plain, plainMap);

        // Same k-mer assignment, but with the additional node 7 inserted between 3 and its children,
        // as a refinement would introduce it.
        SmallTaxTree refined = buildTree(new int[][]{{1, 1}, {2, 1}, {3, 2}, {4, 2}, {7, 3}, {5, 7}, {6, 7}});
        Map<String, DBQualityCountsGoal.Counts> refinedMap = baseCounts();
        refinedMap.put("7", counts("7", 0, 0));
        DBQualityCountsGoal.aggregateCounts(refined, refinedMap);

        DBQualityCountsGoal.Counts plainM = plainMap.get("2");
        DBQualityCountsGoal.Counts refinedM = refinedMap.get("2");
        assertEquals(plainM.getSubtreeTp(), refinedM.getSubtreeTp());
        assertEquals(plainM.getSubtreeTpPlusFp(), refinedM.getSubtreeTpPlusFp());
        assertEquals(plainM.getPrecision(), refinedM.getPrecision(), 1e-12);
        assertEquals(plainM.getSubtreeKmerSum(), refinedM.getSubtreeKmerSum());
        assertEquals(plainM.getSubtreePrecision(), refinedM.getSubtreePrecision(), 1e-12);

        // Under a Laplace correction - one virtual k-mer per node, counted per data taxon under it -
        // the extra node would have raised the refined tree above the plain one although not a single
        // k-mer moved.
        assertTrue(laplaceCorrected(refinedMap, "2", "3", "4", "5", "6", "7")
                > laplaceCorrected(plainMap, "2", "3", "4", "5", "6"));
    }

    /**
     * Computes what the weighted path precision of the subtree formed by the given nodes would be if it
     * were Laplace-corrected by one virtual k-mer per node and data taxon under it.
     */
    private double laplaceCorrected(Map<String, DBQualityCountsGoal.Counts> map, String subtreeRoot, String... others) {
        long virtual = map.get(subtreeRoot).getLeaves();
        for (String taxId : others) {
            virtual += map.get(taxId).getLeaves();
        }
        DBQualityCountsGoal.Counts m = map.get(subtreeRoot);
        return ((double) (virtual + m.getSubtreeTp())) / (virtual + m.getSubtreeTpPlusFp());
    }

    @Test
    public void testSubtreeWithoutKMersLeavesTheMeasuresUndefined() throws IOException {
        SmallTaxTree tree = buildTree(new int[][]{{1, 1}, {2, 1}, {3, 2}, {4, 2}, {5, 3}, {6, 3}});
        Map<String, DBQualityCountsGoal.Counts> map = new HashMap<>();
        for (String taxId : new String[]{"1", "2", "3", "4", "5", "6"}) {
            map.put(taxId, counts(taxId, 0, 0));
        }
        DBQualityCountsGoal.aggregateCounts(tree, map);

        DBQualityCountsGoal.Counts m = map.get("2");
        assertEquals(0, m.getSubtreeTpPlusFp());
        assertEquals(0, m.getSubtreeKmerSum());
        assertTrue(Double.isNaN(m.getPrecision()));
        assertTrue(Double.isNaN(m.getNodePrecision()));
        assertTrue(Double.isNaN(m.getSubtreePrecision()));
    }
}
