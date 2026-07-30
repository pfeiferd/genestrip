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
package org.metagene.genestrip.store;

import org.junit.Test;
import org.metagene.genestrip.tax.SmallTaxTree;
import org.metagene.genestrip.tax.SmallTaxTree.SmallTaxIdNode;
import org.metagene.genestrip.tax.TaxTree;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;

import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.Arrays;
import java.util.BitSet;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertSame;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

/**
 * Tests the properties of {@link SmallTaxIdNode#getStoreIndex()} that callers rely on when they use a
 * node's store index as a bit index instead of the node itself, as {@code KMerIndexBloomGoal} does:
 * after {@link Database#initStoreIndices()} every node carrying a taxid must have a non-negative
 * index, and no two nodes may share one. Unlike the node position, the store index is
 * {@code transient} and therefore only established by that call - an unassigned one reads as 0, not
 * as -1, which would silently alias the node holding index 0.
 */
public class DatabaseStoreIndexTest {
	private static final String[] IDS = { "1", "2", "3", "4", "5", "6", "7" };

	private SmallTaxTree buildSmallTree() throws IOException {
		File dir = Files.createTempDirectory("storeindex").toFile();
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

	/**
	 * Builds a database over the given tree. Only a subset of the tax ids is registered up front, so
	 * that the test also covers the ids {@link Database#initStoreIndices()} has to add itself.
	 */
	private Database buildDatabase(SmallTaxTree tree) {
		KMerStore<String> store = new KMerSortedArray<>(2, 0.0001, 0.0001, Arrays.asList("1", "2"), false, true, 3);
		return new Database(store, tree, null);
	}

	/**
	 * Verifies that every node carrying a taxid ends up with a distinct, non-negative store index -
	 * the precondition {@code KMerIndexBloomGoal} guards on before indexing by store index.
	 */
	@Test
	public void testStoreIndicesUsableAsBitIndex() throws IOException {
		SmallTaxTree tree = buildSmallTree();
		Database db = buildDatabase(tree);
		db.initStoreIndices();

		Map<Integer, String> byIndex = new HashMap<>();
		for (SmallTaxIdNode node : tree) {
			if (node.taxId == null) {
				continue;
			}
			int storeIndex = node.getStoreIndex();
			assertTrue("Store index of " + node.taxId + " must not be negative but is " + storeIndex,
					storeIndex >= 0);
			String clash = byIndex.put(storeIndex, node.taxId);
			if (clash != null) {
				fail("Tax ids " + clash + " and " + node.taxId + " share store index " + storeIndex);
			}
		}
		assertEquals(IDS.length, byIndex.size());
	}

	/**
	 * Verifies that a store index maps straight back to its own node through the converted store,
	 * which is how {@code KMerIndexBloomGoal} checks that an index identifies its node. This relies on
	 * {@link KMerStore#convertValues} preserving the index mapping of the store the indices were
	 * assigned from.
	 */
	@Test
	public void testStoreIndexMapsBackToItsNode() throws IOException {
		SmallTaxTree tree = buildSmallTree();
		Database db = buildDatabase(tree);
		db.initStoreIndices();

		KMerStore<SmallTaxIdNode> converted = db.convertKMerStore();
		int nValues = converted.getNValues();
		for (SmallTaxIdNode node : tree) {
			int storeIndex = node.getStoreIndex();
			assertTrue("Store index of " + node.taxId + " out of range: " + storeIndex,
					storeIndex >= 0 && storeIndex < nValues);
			assertSame("Store index of " + node.taxId + " does not map back to its node", node,
					converted.getValueForIndex(storeIndex));
		}
	}

	/**
	 * Verifies that the store indices stay within the store's value-index range, which is what keeps a
	 * bit set addressed by them small.
	 */
	@Test
	public void testStoreIndicesStayInValueRange() throws IOException {
		SmallTaxTree tree = buildSmallTree();
		Database db = buildDatabase(tree);
		db.initStoreIndices();

		int nValues = db.getKmerStore().getNValues();
		assertTrue("Expected at least one value per tax id", nValues >= IDS.length);
		assertTrue("Value count must stay within the store's limit", nValues <= KMerSortedArray.MAX_VALUES);
		for (SmallTaxIdNode node : tree) {
			if (node.taxId != null) {
				assertTrue("Store index of " + node.taxId + " out of range: " + node.getStoreIndex(),
						node.getStoreIndex() < nValues);
			}
		}
	}

	/**
	 * Verifies that marking a subset of nodes in a store-index-addressed bit set answers membership
	 * exactly like an identity-based set does - the substitution {@code KMerIndexBloomGoal} makes.
	 */
	@Test
	public void testBitSetMatchesIdentitySet() throws IOException {
		SmallTaxTree tree = buildSmallTree();
		Database db = buildDatabase(tree);
		db.initStoreIndices();

		// Arbitrary but non-trivial subset: every node whose tax id is an odd number.
		BitSet bits = new BitSet();
		Set<SmallTaxIdNode> expected = new HashSet<>();
		for (SmallTaxIdNode node : tree) {
			if (Integer.parseInt(node.taxId) % 2 == 1) {
				expected.add(node);
				bits.set(node.getStoreIndex());
			}
		}
		assertTrue(expected.size() > 1 && expected.size() < IDS.length);

		for (SmallTaxIdNode node : tree) {
			assertEquals("Mismatch for tax id " + node.taxId, expected.contains(node),
					bits.get(node.getStoreIndex()));
		}
	}
}
