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
package org.metagene.genestrip.tax;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertSame;

import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;

import org.junit.Test;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;

/**
 * Unit tests for the depth-based {@link TaxTree#getLowestCommonAncestor} (optimization 5b), on a
 * small hand-built taxonomy:
 *
 * <pre>
 *   1                      depth 0
 *   +-- 2                  depth 1
 *   |   +-- 3              depth 2
 *   |   |   +-- 5          depth 3
 *   |   |       +-- 6      depth 4
 *   |   +-- 4              depth 2
 *   +-- 7                  depth 1
 * </pre>
 */
public class TaxTreeLCATest {
	private TaxTree buildTree() throws IOException {
		File dir = Files.createTempDirectory("taxtree").toFile();
		dir.deleteOnExit();
		int[][] edges = { { 1, 1 }, { 2, 1 }, { 3, 2 }, { 4, 2 }, { 5, 3 }, { 6, 5 }, { 7, 1 } };
		StringBuilder nodes = new StringBuilder();
		StringBuilder names = new StringBuilder();
		for (int[] e : edges) {
			nodes.append(e[0]).append("\t|\t").append(e[1]).append("\t|\tno rank\t|\t\t|\n");
			names.append(e[0]).append("\t|\t").append(e[0]).append("\t|\t\t|\tscientific name\t|\n");
		}
		Files.write(new File(dir, TaxTree.NODES_DMP).toPath(), nodes.toString().getBytes(StandardCharsets.UTF_8));
		Files.write(new File(dir, TaxTree.NAMES_DMP).toPath(), names.toString().getBytes(StandardCharsets.UTF_8));
		return new TaxTree(dir, false);
	}

	@Test
	public void testDepthsAssigned() throws IOException {
		TaxTree t = buildTree();
		assertEquals(0, t.getNodeByTaxId("1").getDepth());
		assertEquals(1, t.getNodeByTaxId("2").getDepth());
		assertEquals(2, t.getNodeByTaxId("3").getDepth());
		assertEquals(2, t.getNodeByTaxId("4").getDepth());
		assertEquals(3, t.getNodeByTaxId("5").getDepth());
		assertEquals(4, t.getNodeByTaxId("6").getDepth());
		assertEquals(1, t.getNodeByTaxId("7").getDepth());
	}

	@Test
	public void testLowestCommonAncestor() throws IOException {
		TaxTree t = buildTree();
		// Equal nodes.
		assertSame(node(t, "6"), t.getLowestCommonAncestor(node(t, "6"), node(t, "6")));
		// Ancestor / descendant (either order).
		assertSame(node(t, "5"), t.getLowestCommonAncestor(node(t, "5"), node(t, "6")));
		assertSame(node(t, "5"), t.getLowestCommonAncestor(node(t, "6"), node(t, "5")));
		assertSame(node(t, "1"), t.getLowestCommonAncestor(node(t, "1"), node(t, "6")));
		// Diverging within a subtree: 6 (under 3 under 2) and 4 (under 2) meet at 2.
		assertSame(node(t, "2"), t.getLowestCommonAncestor(node(t, "6"), node(t, "4")));
		assertSame(node(t, "2"), t.getLowestCommonAncestor(node(t, "4"), node(t, "6")));
		assertSame(node(t, "2"), t.getLowestCommonAncestor(node(t, "3"), node(t, "4")));
		// Different top-level branches meet at the root.
		assertSame(node(t, "1"), t.getLowestCommonAncestor(node(t, "6"), node(t, "7")));
		// Null handling.
		assertNull(t.getLowestCommonAncestor(null, node(t, "6")));
		assertNull(t.getLowestCommonAncestor(node(t, "6"), null));
	}

	// Cross-checks getLowestCommonAncestor against a brute-force ancestor-set intersection over all
	// ordered node pairs, so the depth-based walk is verified against an independent oracle.
	@Test
	public void testLowestCommonAncestorMatchesBruteForce() throws IOException {
		TaxTree t = buildTree();
		String[] ids = { "1", "2", "3", "4", "5", "6", "7" };
		for (String x : ids) {
			for (String y : ids) {
				TaxIdNode a = node(t, x);
				TaxIdNode b = node(t, y);
				assertSame("LCA(" + x + "," + y + ")", bruteForceLCA(a, b), t.getLowestCommonAncestor(a, b));
			}
		}
	}

	private static TaxIdNode bruteForceLCA(TaxIdNode a, TaxIdNode b) {
		for (TaxIdNode x = a; x != null; x = x.getParent()) {
			for (TaxIdNode y = b; y != null; y = y.getParent()) {
				if (x == y) {
					return x;
				}
			}
		}
		return null;
	}

	private static TaxIdNode node(TaxTree t, String id) {
		return t.getNodeByTaxId(id);
	}
}
