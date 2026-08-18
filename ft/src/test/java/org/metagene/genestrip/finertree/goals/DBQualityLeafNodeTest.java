package org.metagene.genestrip.finertree.goals;

import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.LinkedHashMap;
import java.util.Map;

import org.junit.Test;
import org.metagene.genestrip.tax.Rank;
import org.metagene.genestrip.tax.SmallTaxTree;
import org.metagene.genestrip.tax.TaxTree;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

/**
 * Tests {@link DBQualityCountsGoal#isLeafNode}, which decides what the intrinsic quality measures are
 * taken over.
 * <p>
 * The rule has to agree with both halves of the goal at once. The database fill nests the artificial
 * nodes tax id -&gt; DATA -&gt; FILE -&gt; ID and stores a genomic file's k-mers at the deepest of them
 * that is enabled; reading the fastas back, the reader resolves a record to the deepest that exists.
 * A leaf is therefore the deepest artificial node on its branch, and testing for {@link Rank#DATA}
 * alone -- which is what this did before -- is right only while file and id nodes are both off.
 */
public class DBQualityLeafNodeTest {

    /**
     * Builds a {@link SmallTaxTree} from {@code {child, parent}} edges, giving each node the rank named
     * for it in {@code ranks} and "no rank" otherwise.
     */
    private SmallTaxTree buildTree(int[][] edges, Map<String, Rank> ranks) throws IOException {
        File dir = Files.createTempDirectory("dbqualityleaf").toFile();
        dir.deleteOnExit();
        StringBuilder nodes = new StringBuilder();
        StringBuilder names = new StringBuilder();
        for (int[] e : edges) {
            Rank rank = ranks.get(String.valueOf(e[0]));
            nodes.append(e[0]).append("\t|\t").append(e[1]).append("\t|\t")
                 .append(rank == null ? "no rank" : rank.getName()).append("\t|\t\t|\n");
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

    private boolean isLeaf(SmallTaxTree tree, String taxId) {
        SmallTaxTree.SmallTaxIdNode node = tree.getNodeByTaxId(taxId);
        assertTrue("node " + taxId + " is missing from the tree", node != null);
        return DBQualityCountsGoal.isLeafNode(node);
    }

    /**
     * The refined shape, and the case that used to be got wrong. `ftupdatedb' inserts REFINED nodes
     * BETWEEN the data node and its file nodes, so the file nodes stop being its children. Judged by
     * the direct children alone the data node then looks like the deepest artificial node on its
     * branch and counts as a leaf -- whereupon its k-mers, the ones a refinement has the most to gain
     * on, drop out of every average restricted to what sits above the data taxa, and the same
     * database reports a different denominator before and after being refined. On `cdiff' that alone
     * lifted the reported sp* from 0.236 to 0.338 with no k-mer moving at all.
     */
    @Test
    public void aDataNodeStaysInternalWhenRefinedNodesSitBelowIt() throws IOException {
        Map<String, Rank> ranks = new LinkedHashMap<>();
        ranks.put("3", Rank.DATA);
        ranks.put("4", Rank.REFINED);
        ranks.put("5", Rank.REFINED);
        ranks.put("6", Rank.FILE);
        ranks.put("7", Rank.FILE);
        // 1 -> 2 -> 3(DATA) -> 4(REFINED) -> 5(REFINED) -> {6(FILE), 7(FILE)}
        SmallTaxTree tree = buildTree(new int[][] { { 1, 1 }, { 2, 1 }, { 3, 2 }, { 4, 3 }, { 5, 4 },
                { 6, 5 }, { 7, 5 } }, ranks);
        assertFalse("the data node is internal: file nodes sit below it, two refined levels down",
                isLeaf(tree, "3"));
        assertFalse("a refined node is never a leaf, it has no origin rank", isLeaf(tree, "4"));
        assertFalse("a refined node is never a leaf, it has no origin rank", isLeaf(tree, "5"));
        assertTrue("the file nodes are the leaves", isLeaf(tree, "6"));
        assertTrue("the file nodes are the leaves", isLeaf(tree, "7"));
    }

    /**
     * A node with children is never a leaf, whatever those children are. The rule asks for the
     * children themselves rather than for their ranks, so nothing the refinement inserts between a
     * node and its original children can make it look like the end of a branch.
     */
    @Test
    public void aDataNodeWithRefinedChildrenIsNotALeaf() throws IOException {
        Map<String, Rank> ranks = new LinkedHashMap<>();
        ranks.put("3", Rank.DATA);
        ranks.put("4", Rank.REFINED);
        ranks.put("5", Rank.REFINED);
        SmallTaxTree tree = buildTree(new int[][] { { 1, 1 }, { 2, 1 }, { 3, 2 }, { 4, 3 }, { 5, 4 } }, ranks);
        assertFalse(isLeaf(tree, "3"));
    }

    /**
     * With file nodes below the data node, the file nodes are the leaves and the data node is an empty
     * intermediate. This is the `cdiff' shape: the taxonomy offers no children below the species, so
     * every assembly becomes a file node and the refinement clusters those.
     */
    @Test
    public void fileNodesAreTheLeavesAndTheirDataParentIsNot() throws IOException {
        Map<String, Rank> ranks = new LinkedHashMap<>();
        ranks.put("3", Rank.DATA);
        ranks.put("4", Rank.FILE);
        ranks.put("5", Rank.FILE);
        //  1 -> 2 -> 3 (DATA) -> {4 (FILE), 5 (FILE)}
        SmallTaxTree tree = buildTree(new int[][]{{1, 1}, {2, 1}, {3, 2}, {4, 3}, {5, 3}}, ranks);

        assertTrue("a file node is where the k-mers come to rest", isLeaf(tree, "4"));
        assertTrue(isLeaf(tree, "5"));
        assertFalse("a data node with file children holds no k-mers of its own", isLeaf(tree, "3"));
        assertFalse("a taxonomy node is never a leaf", isLeaf(tree, "2"));
        assertFalse(isLeaf(tree, "1"));
    }

    /**
     * The behaviour every other database of the paper relies on, unchanged: with no file or id nodes,
     * a data node is the deepest artificial node and therefore the leaf.
     */
    @Test
    public void aDataNodeWithoutFileChildrenIsStillTheLeaf() throws IOException {
        Map<String, Rank> ranks = new LinkedHashMap<>();
        ranks.put("3", Rank.DATA);
        ranks.put("4", Rank.DATA);
        //  1 -> 2 -> {3 (DATA), 4 (DATA)}
        SmallTaxTree tree = buildTree(new int[][]{{1, 1}, {2, 1}, {3, 2}, {4, 2}}, ranks);

        assertTrue(isLeaf(tree, "3"));
        assertTrue(isLeaf(tree, "4"));
        assertFalse(isLeaf(tree, "2"));
    }

    /**
     * An id node below a file node takes the leaf role from it, since that is where the fill would put
     * the k-mers and where the reader would resolve a record to.
     */
    @Test
    public void anIdNodeTakesTheLeafRoleFromItsFileParent() throws IOException {
        Map<String, Rank> ranks = new LinkedHashMap<>();
        ranks.put("3", Rank.DATA);
        ranks.put("4", Rank.FILE);
        ranks.put("5", Rank.ID);
        //  1 -> 2 -> 3 (DATA) -> 4 (FILE) -> 5 (ID)
        SmallTaxTree tree = buildTree(new int[][]{{1, 1}, {2, 1}, {3, 2}, {4, 3}, {5, 4}}, ranks);

        assertTrue(isLeaf(tree, "5"));
        assertFalse(isLeaf(tree, "4"));
        assertFalse(isLeaf(tree, "3"));
    }

    /**
     * A refined node is not an origin rank. The refinement inserts it *above* the artificial nodes and
     * moves k-mers onto it, so it is internal in the way a taxonomy node is -- and the measures
     * restricted to what sits above the data have to keep counting it. Were it treated as a leaf, the
     * k-mers a refinement moves would be excluded from exactly the measure meant to show the movement.
     */
    @Test
    public void aRefinedNodeIsNeverALeaf() throws IOException {
        Map<String, Rank> ranks = new LinkedHashMap<>();
        ranks.put("3", Rank.REFINED);
        ranks.put("4", Rank.DATA);
        //  1 -> 2 -> 3 (REFINED) -> 4 (DATA)
        SmallTaxTree tree = buildTree(new int[][]{{1, 1}, {2, 1}, {3, 2}, {4, 3}}, ranks);

        assertFalse("a refined node holds k-mers but is not where a file's k-mers rest", isLeaf(tree, "3"));
        assertTrue(isLeaf(tree, "4"));
    }
}
