/*
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 *
 * License: Apache 2.0
 *
 * Licensor: Daniel Pfeifer (daniel.pfeifer@progotec.de)
 *
 */
package org.metagene.genestrip.tax;

import org.junit.Before;
import org.junit.Test;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;

import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertSame;
import static org.junit.Assert.assertTrue;

/**
 * Tests {@link TaxTree#otherNode(TaxIdNode, TaxTree.IDStringGenerator)}, the artificial node that
 * stands for the genomes below a taxon which nothing in the database names.
 * <p>
 * The properties under test are the ones the database build relies on: the node is created once and
 * found again by rank, it is a sibling of the taxon's {@code DATA} node rather than a step in the
 * {@code DATA} - {@code FILE} - {@code GENOME} - {@code ID} chain, and a lookup without a generator
 * never creates one - which is what lets the fill pass look up what the counting pass created.
 */
public class OtherNodeTest {
	private TaxTree tree;
	private int counter;

	private TaxTree buildTree() throws IOException {
		File dir = Files.createTempDirectory("othernode").toFile();
		dir.deleteOnExit();
		int[][] edges = { { 1, 1 }, { 2, 1 }, { 3, 1 } };
		StringBuilder nodes = new StringBuilder();
		StringBuilder names = new StringBuilder();
		for (int[] e : edges) {
			nodes.append(e[0]).append("\t|\t").append(e[1]).append("\t|\tspecies\t|\t\t|\n");
			names.append(e[0]).append("\t|\t").append(e[0]).append("\t|\t\t|\tscientific name\t|\n");
		}
		Files.write(new File(dir, TaxTree.NODES_DMP).toPath(), nodes.toString().getBytes(StandardCharsets.UTF_8));
		Files.write(new File(dir, TaxTree.NAMES_DMP).toPath(), names.toString().getBytes(StandardCharsets.UTF_8));
		return new TaxTree(dir, false);
	}

	/** Artificial tax ids are digits starting with {@code 00}: they are kept in a digit trie. */
	private TaxIdNode create(TaxIdNode parent) {
		return tree.otherNode(parent, c -> "00" + (counter++));
	}

	@Before
	public void setUp() throws IOException {
		tree = buildTree();
		counter = 0;
	}

	@Test
	public void createdNodeCarriesTheOtherRank() {
		TaxIdNode other = create(tree.getNodeByTaxId("2"));
		assertNotNull(other);
		assertEquals(Rank.OTHER.ordinal(), other.getRankOrdinal());
	}

	/** The name says which taxon the node belongs to, as the refinement's OTHER slot does. */
	@Test
	public void namedAfterItsTaxon() {
		assertEquals("2OTHER", create(tree.getNodeByTaxId("2")).getName());
	}

	@Test
	public void creationIsIdempotent() {
		TaxIdNode parent = tree.getNodeByTaxId("2");
		assertSame(create(parent), create(parent));
		assertEquals(1, counter);
	}

	@Test
	public void theNodeIsFoundAgainByRank() {
		TaxIdNode parent = tree.getNodeByTaxId("2");
		assertSame(create(parent), parent.getOtherChild());
	}

	@Test
	public void absentBeforeCreation() {
		assertNull(tree.getNodeByTaxId("3").getOtherChild());
	}

	@Test
	public void lookUpWithoutGeneratorCreatesNothing() {
		TaxIdNode parent = tree.getNodeByTaxId("3");
		assertNull(tree.otherNode(parent, null));
		assertNull(parent.getOtherChild());
	}

	@Test
	public void eachTaxonGetsItsOwn() {
		TaxIdNode a = create(tree.getNodeByTaxId("2"));
		TaxIdNode b = create(tree.getNodeByTaxId("3"));
		assertTrue(a != b);
		assertEquals(2, counter);
	}

	/**
	 * The OTHER node hangs off the taxon beside its DATA node, not below it: the two are siblings,
	 * which is what keeps OTHER out of the origin chain.
	 */
	@Test
	public void siblingOfTheDataNode() {
		TaxIdNode parent = tree.getNodeByTaxId("2");
		TaxIdNode data = tree.dataNode(parent, c -> "00" + (counter++));
		TaxIdNode other = create(parent);
		assertSame(data, parent.getDataChild());
		assertSame(other, parent.getOtherChild());
		assertNull(other.getDataChild());
		assertEquals(2, parent.getSubNodes().size());
	}
}
