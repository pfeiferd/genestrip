package org.metagene.genestrip.refseq;

import org.junit.Before;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TemporaryFolder;
import org.metagene.genestrip.refseq.AbstractRefSeqFastaReader.StringLong2DigitTrie;
import org.metagene.genestrip.tax.TaxTree;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;
import org.metagene.genestrip.util.ByteArrayUtil;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

/**
 * What a genome is, and that a reader lets a contig in exactly when its genome is in the selection.
 * <p>
 * Counting fasta entries cannot express "fifty genomes": a draft assembly arrives as one entry per
 * contig, so for a species RefSeq holds as drafts a limit of fifty entries is reached inside a single
 * assembly. What ties the contigs of one assembly together is the WGS accession.
 * <p>
 * The selection itself is not made here. It is made once, while the accession catalog is read, and
 * travels with the accession map - see {@link AccessionMap#getAdmittedGenomes()}. A reader only asks.
 */
public class GenomeLimitTest {
	@Rule
	public final TemporaryFolder folder = new TemporaryFolder();

	private TaxTree tree;

	/**
	 * <pre>
	 * 1 root
	 *  +- 1301 genus Streptococcus
	 *       +- 1313 species S. pneumoniae      &lt;- WGS project AAAA, two assembly versions
	 *       |    +- 170187 strain TIGR4        &lt;- WGS project BBBB
	 *       +- 1309 species S. mutans          &lt;- finished replicon NZ_CP000001
	 * </pre>
	 */
	@Before
	public void setUp() throws IOException {
		File dir = folder.newFolder("taxonomy");
		write(new File(dir, "nodes.dmp"),
				"1\t|\t1\t|\tno rank\t|",
				"1301\t|\t1\t|\tgenus\t|",
				"1313\t|\t1301\t|\tspecies\t|",
				"170187\t|\t1313\t|\tstrain\t|",
				"1309\t|\t1301\t|\tspecies\t|");
		write(new File(dir, "names.dmp"),
				"1\t|\troot\t|\t\t|\tscientific name\t|",
				"1301\t|\tStreptococcus\t|\t\t|\tscientific name\t|",
				"1313\t|\tStreptococcus pneumoniae\t|\t\t|\tscientific name\t|",
				"170187\t|\tStreptococcus pneumoniae TIGR4\t|\t\t|\tscientific name\t|",
				"1309\t|\tStreptococcus mutans\t|\t\t|\tscientific name\t|");
		tree = new TaxTree(dir, false);
	}

	// ---- what a genome key is -------------------------------------------------------------------

	/** All contigs of one WGS project share the letter prefix, which is where the key is cut. */
	@Test
	public void wgsContigsOfOneProjectShareTheirKey() {
		assertEquals("NZ_CABEIU", key("NZ_CABEIU010000001.1"));
		assertEquals("NZ_CABEIU", key("NZ_CABEIU010000273.1"));
		assertEquals("NZ_AAAA", key("NZ_AAAA01000001.1"));
	}

	/**
	 * The assembly version stays out of the key. Counting {@code CABEIU01} and {@code CABEIU02} apart
	 * doubles the tally - measured against RefSeq release 233 it turns 8,919 genomes of
	 * <em>S. pneumoniae</em> into 16,771, against the 9,263 assemblies actually on record.
	 */
	@Test
	public void wgsVersionsAreOneGenome() {
		assertEquals(key("NZ_CABEIU010000001.1"), key("NZ_CABEIU020000001.1"));
	}

	/** The WGS master record of a project belongs to that project. */
	@Test
	public void wgsMasterRecordKeepsItsProject() {
		assertEquals(key("NZ_CABEIU010000001.1"), key("NZ_CABEIU000000000.1"));
	}

	/** A finished replicon names itself - minus its version, which is not a different genome. */
	@Test
	public void finishedRepliconsStandForThemselves() {
		assertEquals("NZ_CP012345", key("NZ_CP012345.1"));
		assertEquals("NZ_CP012345", key("NZ_CP012345.2"));
		assertEquals("NC_003028", key("NC_003028.3"));
		assertEquals("X", key("X"));
	}

	/** Two plasmids of one finished assembly do count as two, which is the error we accept. */
	@Test
	public void plasmidsOfAFinishedGenomeCountSeparately() {
		assertEquals(2, new HashSet<>(Arrays.asList(key("NZ_CP012345.1"), key("NZ_CP012346.1"))).size());
	}

	// ---- what the selection holds ---------------------------------------------------------------

	/** A genome counts once however many accessions it arrives in - that is what the set is for. */
	@Test
	public void oneGenomeIsAdmittedOnce() {
		GenomeKeyTrie selection = new GenomeKeyTrie();
		assertTrue(admit(selection, "NZ_AAAA01000001.1"));
		assertFalse(admit(selection, "NZ_AAAA01000002.1"));
		assertFalse(admit(selection, "NZ_AAAA02000001.1"));
		assertTrue(admit(selection, "NZ_BBBB01000001.1"));
	}

	/**
	 * Once the catalog scan has closed the selection, nothing may be added to it: it is read by every
	 * pass and every reader thread from then on, and an addition would be a race and a silent change of
	 * which genomes the database holds.
	 */
	@Test
	public void aClosedSelectionRefusesFurtherGenomes() {
		GenomeKeyTrie selection = selectionOf("NZ_AAAA01000001.1");
		selection.freeze();
		assertFalse(admit(selection, "NZ_AAAA01000002.1"));
		try {
			admit(selection, "NZ_BBBB01000001.1");
			fail("expected the closed selection to refuse a new genome");
		}
		catch (IllegalStateException expected) {
			// The point of the test.
		}
	}

	// ---- what a reader does with it -------------------------------------------------------------

	/** A genome in the selection comes in whole, whatever order its contigs arrive in. */
	@Test
	public void anAdmittedGenomeComesInWhole() throws IOException {
		assertEquals(Arrays.asList("NZ_AAAA01000001.1", "NZ_AAAA01000002.1", "NZ_AAAA02000001.1"),
				read(selectionOf("NZ_AAAA01000001.1")).included);
	}

	/** And a genome that is not in it stays out whole. */
	@Test
	public void aGenomeOutsideTheSelectionStaysOut() throws IOException {
		assertEquals(Arrays.asList("NZ_BBBB01000001.1", "NZ_CP000001.1"),
				read(selectionOf("NZ_BBBB01000001.1", "NZ_CP000001.1")).included);
	}

	/** An empty selection admits nothing at all. */
	@Test
	public void anEmptySelectionAdmitsNothing() throws IOException {
		assertEquals(0, read(new GenomeKeyTrie()).included.size());
	}

	/**
	 * Without a selection the reader is handed none - that is what the accession map returns where
	 * {@code maxGenomesPerTaxid} is not set - and reads everything.
	 */
	@Test
	public void withoutASelectionEverythingIsRead() throws IOException {
		Reader reader = read(null);
		assertEquals(5, reader.included.size());
		assertEquals(5, contigs(reader, "1301"));
	}

	/**
	 * And a map built for a database that is to hold every genome carries no selection at all - there
	 * is no empty set standing in for "no limit", so nothing of the machinery exists to be reached.
	 */
	@Test
	public void aMapWithoutALimitCarriesNoSelection() {
		assertNull(new AccessionMapTrieImpl().getAdmittedGenomes());
		assertNull(new AccessionMapImpl().getAdmittedGenomes());
	}

	// ---- harness --------------------------------------------------------------------------------

	private static String key(String accession) {
		byte[] bytes = accession.getBytes(StandardCharsets.US_ASCII);
		return new String(bytes, 0, GenomeKeyTrie.genomeKeyLength(bytes, 0, bytes.length),
				StandardCharsets.US_ASCII);
	}

	private static boolean admit(GenomeKeyTrie selection, String accession) {
		byte[] bytes = accession.getBytes(StandardCharsets.US_ASCII);
		return selection.admit(bytes, 0, bytes.length);
	}

	private static GenomeKeyTrie selectionOf(String... accessions) {
		GenomeKeyTrie selection = new GenomeKeyTrie();
		for (String accession : accessions) {
			admit(selection, accession);
		}
		return selection;
	}

	private static long contigs(Reader reader, String taxId) {
		return reader.getContigsPerTaxid().get(taxId).getLongValue();
	}

	/**
	 * Reads a fasta of five contigs - three of one WGS project in two assembly versions, one of a
	 * second project filed at a strain below the same species, and one finished replicon of another
	 * species - against the given selection.
	 */
	private Reader read(GenomeKeyTrie selection) throws IOException {
		File fasta = folder.newFile();
		write(fasta,
				">NZ_AAAA01000001.1 first project, first contig", "ACGTACGTACGT",
				">NZ_AAAA01000002.1 first project, second contig", "ACGTACGTACGT",
				">NZ_AAAA02000001.1 first project, next assembly version", "ACGTACGTACGT",
				">NZ_BBBB01000001.1 second project, filed at the strain", "ACGTACGTACGT",
				">NZ_CP000001.1 a finished replicon of another species", "ACGTACGTACGT");
		Reader reader = new Reader(tree, selection);
		reader.readFasta(fasta);
		return reader;
	}

	private static void write(File file, String... lines) throws IOException {
		try (PrintWriter pw = new PrintWriter(file, StandardCharsets.UTF_8.name())) {
			for (String line : lines) {
				pw.println(line);
			}
		}
	}

	/** Resolves an accession by its letters and carries the selection, as the real map does. */
	private static final class LetterAccessionMap implements AccessionMap {
		private final TaxTree tree;
		private final GenomeKeyTrie admittedGenomes;

		private LetterAccessionMap(TaxTree tree, GenomeKeyTrie admittedGenomes) {
			this.tree = tree;
			this.admittedGenomes = admittedGenomes;
		}

		@Override
		public void put(byte[] array, int start, int end, TaxIdNode node) {
		}

		@Override
		public TaxIdNode get(byte[] array, int start, int end, boolean assemblyAccessionsOnly) {
			String accession = new String(array, start, end - start, StandardCharsets.US_ASCII);
			if (accession.startsWith("NZ_AAAA")) {
				return tree.getNodeByTaxId("1313");
			}
			if (accession.startsWith("NZ_BBBB")) {
				return tree.getNodeByTaxId("170187");
			}
			return tree.getNodeByTaxId("1309");
		}

		@Override
		public void optimize() {
		}

		@Override
		public GenomeKeyTrie getAdmittedGenomes() {
			return admittedGenomes;
		}

		@Override
		public int getEntriesForNode(TaxIdNode node) {
			return 0;
		}
	}

	/** Records which contigs were let in and counts one k-mer per data line. */
	private static final class Reader extends AbstractRefSeqFastaReader {
		private final List<String> included = new ArrayList<>();

		private Reader(TaxTree tree, GenomeKeyTrie selection) {
			super(4096, allNodes(tree), new LetterAccessionMap(tree, selection), 31, 1, false,
					new StringLong2DigitTrie());
		}

		private static Set<TaxIdNode> allNodes(TaxTree tree) {
			Set<TaxIdNode> nodes = new HashSet<>();
			for (String taxId : new String[] { "1313", "170187", "1309" }) {
				nodes.add(tree.getNodeByTaxId(taxId));
			}
			return nodes;
		}

		@Override
		protected void infoLine() {
			super.infoLine();
			if (includeContig) {
				int pos = ByteArrayUtil.indexOf(target, 0, size, ' ');
				included.add(new String(target, 1, pos - 1, StandardCharsets.US_ASCII));
			}
		}

		@Override
		protected void dataLine() {
			kmersInContig++;
		}
	}
}
