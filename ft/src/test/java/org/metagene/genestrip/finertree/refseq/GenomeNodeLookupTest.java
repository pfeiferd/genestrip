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
package org.metagene.genestrip.finertree.refseq;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.charset.StandardCharsets;
import java.util.Collections;
import java.util.Set;

import org.junit.Before;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TemporaryFolder;
import org.metagene.genestrip.refseq.AbstractRefSeqFastaReader.StringLong2DigitTrie;
import org.metagene.genestrip.refseq.AccessionMap;
import org.metagene.genestrip.refseq.GenomeKeyTrie;
import org.metagene.genestrip.refseq.ReworkingStoreFastaReader;
import org.metagene.genestrip.tax.Rank;
import org.metagene.genestrip.tax.SmallTaxTree;
import org.metagene.genestrip.tax.TaxTree;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;

/**
 * Tests that {@link AbstractUpdateFastaReader} finds the genome node {@link ReworkingStoreFastaReader}
 * created for a contig, both where the accession names the genome and where the taxon does - see
 * {@code RefSeqCategory#isOneGenomePerTaxon()}.
 */
public class GenomeNodeLookupTest {
	private static final String CHROMOSOME_1 = "NC_004325";
	private static final String CHROMOSOME_2 = "NC_004326";

	@Rule
	public final TemporaryFolder folder = new TemporaryFolder();

	private TaxTree tree;
	private TaxIdNode taxon;

	@Before
	public void setUp() throws IOException {
		File dir = folder.newFolder("taxonomy");
		write(new File(dir, TaxTree.NODES_DMP), "1\t|\t1\t|\tno rank\t|", "2\t|\t1\t|\tspecies\t|");
		write(new File(dir, TaxTree.NAMES_DMP), "1\t|\troot\t|\t\t|\tscientific name\t|",
				"2\t|\tsome species\t|\t\t|\tscientific name\t|");
		tree = new TaxTree(dir, false);
		taxon = tree.getNodeByTaxId("2");
	}

	@Test
	public void testTaxonKeyedGenomeIsFound() throws IOException {
		File file = fasta("protozoa.1.1.genomic.fna");
		new Creator().read(file);
		LookUp lookUp = new LookUp(tree.toSmallTaxTree());
		lookUp.read(file);
		assertNotNull(lookUp.found);
		assertEquals(Rank.GENOME.ordinal(), lookUp.found.getRankOrdinal());
		assertEquals("2", lookUp.found.getName());
	}

	@Test
	public void testTaxonKeyedGenomeIsFoundInComplete() throws IOException {
		taxon.markOneGenomePerTaxon();
		File file = fasta("complete.1.1.genomic.fna");
		new Creator().read(file);
		LookUp lookUp = new LookUp(tree.toSmallTaxTree());
		lookUp.read(file);
		assertNotNull(lookUp.found);
		assertEquals("2", lookUp.found.getName());
	}

	@Test
	public void testAccessionKeyedGenomeIsFound() throws IOException {
		File file = fasta("bacteria.1.1.genomic.fna");
		new Creator().read(file);
		LookUp lookUp = new LookUp(tree.toSmallTaxTree());
		lookUp.read(file);
		assertNotNull(lookUp.found);
		assertEquals(Rank.GENOME.ordinal(), lookUp.found.getRankOrdinal());
		// The last contig read is the second chromosome, which is a genome of its own here.
		assertEquals(CHROMOSOME_2, lookUp.found.getName());
	}

	private File fasta(String name) throws IOException {
		File file = new File(folder.getRoot(), name);
		write(file, ">" + CHROMOSOME_1 + ".1 chromosome 1", "ACGTACGTACGTACGTACGTACGTACGTACGTACGT",
				">" + CHROMOSOME_2 + ".1 chromosome 2", "ACGTACGTACGTACGTACGTACGTACGTACGTACGT");
		return file;
	}

	private static void write(File file, String... lines) throws IOException {
		try (PrintWriter pw = new PrintWriter(file, StandardCharsets.UTF_8.name())) {
			for (String line : lines) {
				pw.println(line);
			}
		}
	}

	/** Creates the genome nodes, as the sizing pass does. */
	private final class Creator extends ReworkingStoreFastaReader {
		private Creator() {
			super(4096, Collections.singleton(taxon), new SingleTaxonMap(), 31, -1, 1, false,
					new StringLong2DigitTrie(), true, tree, false, false, true, false, true,
					counter -> "00" + counter, null);
		}

		private void read(File file) throws IOException {
			ignoreAccessionMap(null);
			readFasta(file);
		}

		@Override
		protected boolean handleStore(long kmer) {
			return false;
		}
	}

	/** Looks the genome nodes up in the small tree, as the ft goals do. */
	private final class LookUp extends AbstractUpdateFastaReader {
		private final SmallTaxTree smallTree;
		private SmallTaxTree.SmallTaxIdNode found;

		private LookUp(SmallTaxTree smallTree) {
			super(4096, Collections.singleton(taxon), new SingleTaxonMap(), 31, -1, 1, false,
					new StringLong2DigitTrie(), true, false, true, false, false);
			this.smallTree = smallTree;
		}

		private void read(File file) throws IOException {
			ignoreAccessionMap(null);
			readFasta(file);
		}

		@Override
		protected SmallTaxTree getTree() {
			return smallTree;
		}

		@Override
		protected void infoLine() {
			super.infoLine();
			found = leafNode;
		}

		@Override
		protected boolean handleStore(long kmer) {
			return false;
		}
	}

	private final class SingleTaxonMap implements AccessionMap {
		@Override
		public void put(byte[] array, int start, int end, TaxIdNode node) {
		}

		@Override
		public TaxIdNode get(byte[] array, int start, int end, boolean assemblyAccessionsOnly) {
			return taxon;
		}

		@Override
		public void optimize() {
		}

		@Override
		public GenomeKeyTrie getAdmittedGenomes() {
			return null;
		}

		@Override
		public int getEntriesForNode(TaxIdNode node) {
			return 0;
		}
	}
}
