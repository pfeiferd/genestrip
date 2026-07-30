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

import org.junit.Test;
import org.metagene.genestrip.tax.SmallTaxTree.SmallTaxIdNode;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.BitSet;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

/**
 * Tests the properties of {@link SmallTaxIdNode#getPosition()} that callers rely on when they use a
 * node's position as an array or bit-vector index instead of the node itself, as
 * {@code KMerIndexBloomGoal} does: positions must be non-negative and unique within a tree, and they
 * must survive serialization, because a tree obtained from a database has been deserialized.
 */
public class SmallTaxTreePositionTest {
	private static final String[] IDS = { "1", "2", "3", "4", "5", "6", "7" };

	private SmallTaxTree buildSmallTree() throws IOException {
		File dir = Files.createTempDirectory("smalltaxtreepos").toFile();
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
		// Mirrors FillDBGoal, which re-establishes the positions of the large tree before deriving the
		// small one - the small tree inherits them rather than renumbering.
		full.reinitPositions();
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

	private void checkPositionsUsableAsIndex(SmallTaxTree tree) {
		Map<Integer, String> byPosition = new HashMap<>();
		int count = 0;
		for (SmallTaxIdNode node : tree) {
			int position = node.getPosition();
			assertTrue("Position of " + node.taxId + " must not be negative but is " + position, position >= 0);
			String clash = byPosition.put(position, node.taxId);
			if (clash != null) {
				fail("Tax ids " + clash + " and " + node.taxId + " share position " + position);
			}
			count++;
		}
		assertEquals(IDS.length, count);
		assertEquals(count, byPosition.size());
	}

	/**
	 * Verifies that the positions of a freshly derived small tree can be used as indices.
	 */
	@Test
	public void testPositionsUsableAsIndex() throws IOException {
		checkPositionsUsableAsIndex(buildSmallTree());
	}

	/**
	 * Verifies that the positions still identify the nodes after a serialization round trip, which is
	 * how a tree reaches the goals that index by position.
	 */
	@Test
	public void testPositionsSurviveSerialization() throws IOException, ClassNotFoundException {
		SmallTaxTree original = buildSmallTree();
		Map<String, Integer> before = new HashMap<>();
		for (SmallTaxIdNode node : original) {
			before.put(node.taxId, node.getPosition());
		}

		SmallTaxTree loaded = roundTrip(original);
		checkPositionsUsableAsIndex(loaded);
		for (SmallTaxIdNode node : loaded) {
			assertEquals("Position of " + node.taxId + " changed", before.get(node.taxId),
					Integer.valueOf(node.getPosition()));
		}
	}

	/**
	 * Verifies that marking a subset of nodes in a position-indexed bit vector answers membership
	 * exactly like an identity-based set does - the substitution {@code KMerIndexBloomGoal} makes.
	 */
	@Test
	public void testBitVectorMatchesIdentitySet() throws IOException {
		SmallTaxTree tree = buildSmallTree();

		// Arbitrary but non-trivial subset: every node whose tax id is an odd number.
		BitSet bits = new BitSet();
		Set<SmallTaxIdNode> expected = new HashSet<>();
		for (SmallTaxIdNode node : tree) {
			if (Integer.parseInt(node.taxId) % 2 == 1) {
				expected.add(node);
				bits.set(node.getPosition());
			}
		}
		assertTrue(expected.size() > 1 && expected.size() < IDS.length);

		for (SmallTaxIdNode node : tree) {
			assertEquals("Mismatch for tax id " + node.taxId, expected.contains(node),
					bits.get(node.getPosition()));
		}
	}
}
