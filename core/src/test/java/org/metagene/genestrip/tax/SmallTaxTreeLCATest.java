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

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;

import org.junit.Test;
import org.metagene.genestrip.tax.SmallTaxTree.SmallTaxIdNode;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;

/**
 * Tests the depth-based {@link SmallTaxTree#getLowestCommonAncestor} (optimization 5b ported to the
 * compact tree), including a serialization round-trip that verifies the transient {@code depth} is
 * re-initialized on load (so the on-disk format is unchanged and old databases keep working).
 *
 * Tree (same shape as {@link TaxTreeLCATest}):
 * <pre>
 *   1 (0) -- 2 (1) -- 3 (2) -- 5 (3) -- 6 (4)
 *   |        |     -- 4 (2)
 *   |     -- 7 (1)
 * </pre>
 */
public class SmallTaxTreeLCATest {
	private static final String[] IDS = { "1", "2", "3", "4", "5", "6", "7" };

	private SmallTaxTree buildSmallTree() throws IOException {
		File dir = Files.createTempDirectory("smalltaxtree").toFile();
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
		TaxTree full = new TaxTree(dir, false);
		for (String id : IDS) {
			TaxIdNode n = full.getNodeByTaxId(id);
			n.markRequired();
		}
		return full.toSmallTaxTree();
	}

	private static SmallTaxTree roundTrip(SmallTaxTree tree) throws IOException, ClassNotFoundException {
		ByteArrayOutputStream bos = new ByteArrayOutputStream();
		try (ObjectOutputStream oos = new ObjectOutputStream(bos)) {
			oos.writeObject(tree);
		}
		try (ObjectInputStream ois = new ObjectInputStream(new ByteArrayInputStream(bos.toByteArray()))) {
			return (SmallTaxTree) ois.readObject();
		}
	}

	private void checkLCAs(SmallTaxTree t) {
		assertSame(node(t, "5"), t.getLowestCommonAncestor(node(t, "5"), node(t, "6")));
		assertSame(node(t, "5"), t.getLowestCommonAncestor(node(t, "6"), node(t, "5")));
		assertSame(node(t, "1"), t.getLowestCommonAncestor(node(t, "1"), node(t, "6")));
		assertSame(node(t, "2"), t.getLowestCommonAncestor(node(t, "6"), node(t, "4")));
		assertSame(node(t, "2"), t.getLowestCommonAncestor(node(t, "3"), node(t, "4")));
		assertSame(node(t, "1"), t.getLowestCommonAncestor(node(t, "6"), node(t, "7")));
		assertSame(node(t, "6"), t.getLowestCommonAncestor(node(t, "6"), node(t, "6")));
		assertNull(t.getLowestCommonAncestor(null, node(t, "6")));
		assertNull(t.getLowestCommonAncestor(node(t, "6"), null));
		// Cross-check against a brute-force ancestor-set intersection over all ordered pairs.
		for (String x : IDS) {
			for (String y : IDS) {
				assertSame("LCA(" + x + "," + y + ")", bruteForceLCA(node(t, x), node(t, y)),
						t.getLowestCommonAncestor(node(t, x), node(t, y)));
			}
		}
	}

	@Test
	public void testLowestCommonAncestorOnBuiltTree() throws Exception {
		checkLCAs(buildSmallTree());
	}

	@Test
	public void testLowestCommonAncestorAfterDeserialization() throws Exception {
		// The whole point: depth is transient, so it must be re-established by readObject on load.
		SmallTaxTree loaded = roundTrip(buildSmallTree());
		checkLCAs(loaded);
	}

	@Test
	public void testGetLevelMatchesDepthAfterLoad() throws Exception {
		SmallTaxTree loaded = roundTrip(buildSmallTree());
		int[] expected = { 0, 1, 2, 2, 3, 4, 1 };
		for (int i = 0; i < IDS.length; i++) {
			assertEquals("level of " + IDS[i], expected[i], node(loaded, IDS[i]).getLevel());
		}
	}

	private static SmallTaxIdNode bruteForceLCA(SmallTaxIdNode a, SmallTaxIdNode b) {
		for (SmallTaxIdNode x = a; x != null; x = x.getParent()) {
			for (SmallTaxIdNode y = b; y != null; y = y.getParent()) {
				if (x == y) {
					return x;
				}
			}
		}
		return null;
	}

	private static SmallTaxIdNode node(SmallTaxTree t, String id) {
		return t.getNodeByTaxId(id);
	}
}
