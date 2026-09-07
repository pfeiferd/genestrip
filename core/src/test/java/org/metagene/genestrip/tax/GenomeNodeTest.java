/*
 * Genestrip
 *
 * Licensor: Daniel Pfeifer (daniel.pfeifer@progotec.de)
 */
package org.metagene.genestrip.tax;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertNotSame;
import static org.junit.Assert.assertSame;

import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;

import org.junit.Test;
import org.metagene.genestrip.refseq.GenomeKeyTrie;
import org.metagene.genestrip.tax.TaxTree.IDStringGenerator;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;

/**
 * Tests {@link TaxTree#genomeNode(TaxIdNode, byte[], int, int, IDStringGenerator)}, the artificial
 * node {@code genomeNodes=true} gives each genome.
 * <p>
 * The property under test is the one the database build depends on: the counting pass creates these
 * nodes and the fill pass only looks them up, by the same key, so a disagreement between the two
 * would file k-mers at a node the store never registered. The key is the genome key of
 * {@link GenomeKeyTrie}, so every contig of one draft assembly must reach the same node and two
 * assemblies must not share one.
 */
public class GenomeNodeTest {
	/** Two contigs of one WGS assembly: same letter prefix, different contig numbers. */
	private static final String CONTIG_A1 = "NZ_CABEIU010000001";
	private static final String CONTIG_A2 = "NZ_CABEIU010000042";
	/** A contig of a different WGS assembly. */
	private static final String CONTIG_B1 = "NZ_CABEIV010000001";
	/** Two finished replicons, which carry no project prefix and stand for themselves. */
	private static final String REPLICON_C = "NC_003098";
	private static final String REPLICON_D = "NC_012469";

	private TaxTree tree;
	private int counter;

	private TaxTree buildTree() throws IOException {
		File dir = Files.createTempDirectory("genomenode").toFile();
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

	/**
	 * Creates the genome node for an accession, as the counting pass does.
	 * <p>
	 * The artificial tax ids are digits and start with {@code 00}, as the generators in the fill
	 * goals produce them: they are kept in a digit trie, which takes nothing else.
	 */
	private TaxIdNode create(TaxIdNode parent, String accession) {
		byte[] b = accession.getBytes(StandardCharsets.UTF_8);
		return tree.genomeNode(parent, b, 0, b.length, c -> "00" + (counter++));
	}

	/** Looks the genome node up by its key, as the fill pass does. */
	private TaxIdNode lookUp(TaxIdNode parent, String accession) {
		byte[] b = accession.getBytes(StandardCharsets.UTF_8);
		return parent.getChildWithName(b, 0, GenomeKeyTrie.genomeKeyLength(b, 0, b.length));
	}

	@Test
	public void contigsOfOneAssemblyShareANode() throws IOException {
		tree = buildTree();
		TaxIdNode taxon = tree.getNodeByTaxId("2");
		assertSame(create(taxon, CONTIG_A1), create(taxon, CONTIG_A2));
	}

	@Test
	public void differentAssembliesGetDifferentNodes() throws IOException {
		tree = buildTree();
		TaxIdNode taxon = tree.getNodeByTaxId("2");
		assertNotSame(create(taxon, CONTIG_A1), create(taxon, CONTIG_B1));
	}

	@Test
	public void finishedRepliconsStandForThemselves() throws IOException {
		tree = buildTree();
		TaxIdNode taxon = tree.getNodeByTaxId("2");
		assertNotSame(create(taxon, REPLICON_C), create(taxon, REPLICON_D));
	}

	@Test
	public void theNodeCarriesTheGenomeRank() throws IOException {
		tree = buildTree();
		TaxIdNode g = create(tree.getNodeByTaxId("2"), CONTIG_A1);
		assertEquals(Rank.GENOME.ordinal(), g.getRankOrdinal());
	}

	/**
	 * The counting pass creates, the fill pass looks up, and the two must land on the same node --
	 * for every contig of the assembly and not only for the one that created it.
	 */
	@Test
	public void fillFindsWhatTheCountingPassCreated() throws IOException {
		tree = buildTree();
		TaxIdNode taxon = tree.getNodeByTaxId("2");
		TaxIdNode created = create(taxon, CONTIG_A1);
		assertNotNull(lookUp(taxon, CONTIG_A1));
		assertSame(created, lookUp(taxon, CONTIG_A1));
		assertSame(created, lookUp(taxon, CONTIG_A2));
		assertSame(create(taxon, REPLICON_C), lookUp(taxon, REPLICON_C));
	}

	/** The node is a child of its taxon, so one assembly key under two taxa is two nodes. */
	@Test
	public void nodesAreKeyedPerTaxon() throws IOException {
		tree = buildTree();
		TaxIdNode one = create(tree.getNodeByTaxId("2"), CONTIG_A1);
		TaxIdNode two = create(tree.getNodeByTaxId("3"), CONTIG_A1);
		assertNotSame(one, two);
		assertSame(one, lookUp(tree.getNodeByTaxId("2"), CONTIG_A1));
		assertSame(two, lookUp(tree.getNodeByTaxId("3"), CONTIG_A1));
	}

	/**
	 * A taxon that filed a genome has a child, which is the property the ft quality goals rest on:
	 * {@code AbstractDBQualityGoal.isLeafNode} calls a childless node a leaf, so no taxonomy node
	 * may be left holding a genome's k-mers with nothing beneath it.
	 */
	@Test
	public void filingAGenomeGivesTheTaxonAChild() throws IOException {
		tree = buildTree();
		TaxIdNode taxon = tree.getNodeByTaxId("2");
		assertNull(taxon.getSubNodes());
		create(taxon, CONTIG_A1);
		assertEquals(1, taxon.getSubNodes().size());
	}
}
