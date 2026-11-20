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
package org.metagene.genestrip.refseq;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.charset.StandardCharsets;
import java.util.Collections;

import org.junit.Before;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TemporaryFolder;
import org.metagene.genestrip.refseq.AbstractRefSeqFastaReader.StringLong2DigitTrie;
import org.metagene.genestrip.tax.Rank;
import org.metagene.genestrip.tax.TaxTree;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;

/**
 * Tests that {@code genomeNodes=true} takes all sequences of one taxon in a eukaryotic release file for
 * one genome (see {@link RefSeqCategory#isOneGenomePerTaxon()}), and tells genomes apart by accession
 * in every other file.
 */
public class EukaryoteGenomeNodesTest {
	/** Two chromosomes of one organism, with consecutive accessions and no shared prefix. */
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
	public void testCategoryFlags() {
		for (RefSeqCategory cat : RefSeqCategory.values()) {
			boolean eukaryote = cat == RefSeqCategory.FUNGI || cat == RefSeqCategory.INVERTEBRATE
					|| cat == RefSeqCategory.PLANT || cat == RefSeqCategory.PROTOZOA
					|| cat == RefSeqCategory.VERTEBRATE_MAMMALIAN || cat == RefSeqCategory.VERTEBRATE_OTHER;
			assertEquals(cat.toString(), eukaryote, cat.isOneGenomePerTaxon());
		}
		assertFalse(RefSeqCategory.COMPLETE.isGenomeKeyed());
		assertEquals(RefSeqCategory.VERTEBRATE_MAMMALIAN,
				RefSeqCategory.fromFileName("vertebrate_mammalian.1.1.genomic.fna.gz"));
	}

	@Test
	public void testProtozoaReleaseFileGetsOneGenomePerTaxon() throws IOException {
		read("protozoa.1.1.genomic.fna", null);
		assertTaxonKeyed();
	}

	@Test
	public void testVertebrateReleaseFileGetsOneGenomePerTaxon() throws IOException {
		read("vertebrate_other.12.1.genomic.fna", null);
		assertTaxonKeyed();
	}

	private void assertTaxonKeyed() {
		TaxIdNode genome = taxon.getChildWithName("2");
		assertNotNull(genome);
		assertEquals(Rank.GENOME.ordinal(), genome.getRankOrdinal());
		assertNull(childNamed(CHROMOSOME_1));
		assertNull(childNamed(CHROMOSOME_2));
	}

	@Test
	public void testOtherReleaseFileGetsGenomePerAccession() throws IOException {
		read("bacteria.1.1.genomic.fna", null);
		assertAccessionKeyed();
	}

	@Test
	public void testCompleteReleaseFileGetsGenomePerAccession() throws IOException {
		read("complete.1.1.genomic.fna", null);
		assertAccessionKeyed();
	}

	@Test
	public void testCompleteReleaseFileGetsOneGenomeForMarkedTaxon() throws IOException {
		taxon.markOneGenomePerTaxon();
		read("complete.1.1.genomic.fna", null);
		assertTaxonKeyed();
	}

	@Test
	public void testMarkOnlyCountsInCompleteFiles() throws IOException {
		// A taxon marked by the catalog still has its genomes told apart in a category's own file.
		taxon.markOneGenomePerTaxon();
		read("mitochondrion.1.1.genomic.fna", null);
		assertAccessionKeyed();
	}

	@Test
	public void testCatalogRows() {
		assertTrue(row("complete|protozoa"));
		assertTrue(row("complete|mitochondrion|vertebrate_other"));
		assertTrue(row("complete|plant|plastid"));
		assertFalse(row("bacteria|complete"));
		assertFalse(row("bacteria|complete|plasmid"));
		assertFalse(row("complete|viral"));
		assertFalse(row("complete|other"));
	}

	private static boolean row(String dirs) {
		byte[] b = ("5833\tPlasmodium\tNC_004325.2\t" + dirs + "\tREVIEWED\t640851")
				.getBytes(StandardCharsets.US_ASCII);
		int start = "5833\tPlasmodium\tNC_004325.2\t".length();
		return RefSeqCategory.isOneGenomePerTaxonRow(b, start, start + dirs.length());
	}

	@Test
	public void testProjectFastaGetsGenomePerAccession() throws IOException {
		// A project's own fasta is not a release file, whatever it is called.
		read("protozoa.1.1.genomic.fna", taxon);
		assertAccessionKeyed();
	}

	@Test
	public void testEukaryoteFileDoesNotStickToTheReader() throws IOException {
		Reader reader = new Reader();
		reader.ignoreAccessionMap(null);
		reader.readFasta(fasta("protozoa.1.1.genomic.fna"));
		reader.readFasta(fasta("bacteria.1.1.genomic.fna"));
		assertNotNull(childNamed(CHROMOSOME_1));
	}

	private void assertAccessionKeyed() {
		assertNotNull(childNamed(CHROMOSOME_1));
		assertNotNull(childNamed(CHROMOSOME_2));
		assertNull(taxon.getChildWithName("2"));
	}

	private TaxIdNode childNamed(String key) {
		byte[] b = key.getBytes(StandardCharsets.US_ASCII);
		return taxon.getChildWithName(b, 0, b.length);
	}

	private void read(String name, TaxIdNode node) throws IOException {
		Reader reader = new Reader();
		reader.ignoreAccessionMap(node);
		reader.readFasta(fasta(name));
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

	/** Creates the genome nodes, as the sizing pass does; resolves every accession to the taxon. */
	private final class Reader extends ReworkingStoreFastaReader {
		private Reader() {
			super(4096, Collections.singleton(taxon), new SingleTaxonMap(), 31, -1, 1, false,
					new StringLong2DigitTrie(), true, tree, false, false, true, false, true,
					counter -> "00" + counter, null);
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
