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
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import static org.junit.Assert.assertEquals;

/**
 * The genome selection over real RefSeq records, in the order the release files hold them.
 * <p>
 * {@code src/test/resources/refseq/bacteria_wgs.fna} is fourteen consecutive records of
 * {@code bacteria.1.1.genomic.fna.gz} (release 233) and {@code refseq/viral_segments.fna} is the
 * eight records of one segmented virus from {@code viral.1.1.genomic.fna.gz} (release 236); only the
 * sequence is truncated, to seventy bases. The tax ids are the ones the RefSeq catalog files those
 * accessions under: 1428 and 1423 for the two bacterial assemblies, 1980428 for all eight segments.
 * <p>
 * The bacterial excerpt is here for one reason: the contigs of an assembly are <em>not</em> gathered
 * in the release. Its fourteen records are two assemblies in five runs, and over the whole file 204
 * assemblies arrive in 4,832 runs. Anything that counted runs rather than assemblies would count each
 * of them two dozen times, and no synthetic fixture would have shown it.
 */
public class GenomeLimitRefSeqTest {
	@Rule
	public final TemporaryFolder folder = new TemporaryFolder();

	private TaxTree tree;

	/** The three taxa the excerpts resolve to, with the tax ids, names and ranks NCBI gives them. */
	@Before
	public void setUp() throws IOException {
		File dir = folder.newFolder("taxonomy");
		write(new File(dir, "nodes.dmp"),
				"1\t|\t1\t|\tno rank\t|",
				"1386\t|\t1\t|\tgenus\t|",
				"1423\t|\t1386\t|\tspecies\t|",
				"1428\t|\t1386\t|\tspecies\t|",
				"1980428\t|\t1\t|\tspecies\t|");
		write(new File(dir, "names.dmp"),
				"1\t|\troot\t|\t\t|\tscientific name\t|",
				"1386\t|\tBacillus\t|\t\t|\tscientific name\t|",
				"1423\t|\tBacillus subtilis\t|\t\t|\tscientific name\t|",
				"1428\t|\tBacillus thuringiensis\t|\t\t|\tscientific name\t|",
				"1980428\t|\tHigh Plains wheat mosaic emaravirus\t|\t\t|\tscientific name\t|");
		tree = new TaxTree(dir, false);
	}

	/**
	 * Twelve contigs of one assembly and two of another, arriving in five runs, are two genomes.
	 * Counting runs would make it five.
	 */
	@Test
	public void interleavedContigsAreOneGenomeEachAssembly() throws IOException {
		assertEquals(2, admitAll("bacteria_wgs.fna"));
	}

	/**
	 * A selection holding only the first assembly keeps every one of its contigs, including the seven
	 * that arrive after the contigs of the other one.
	 */
	@Test
	public void theSelectedAssemblyKeepsItsLaterContigs() throws IOException {
		Reader reader = read("bacteria_wgs.fna", selectionOf("NZ_JARSVH010000408.1"));
		assertEquals(12, reader.included.size());
		for (String accession : reader.included) {
			assertEquals("NZ_JARSVH", accession.substring(0, 9));
		}
	}

	/**
	 * A segmented virus is where the accession stops being evidence of an assembly: the eight segments
	 * of High Plains wheat mosaic emaravirus are eight complete records under one tax id, so they count
	 * as eight genomes. This is the over-count the finished-replicon rule accepts, and it is worth
	 * seeing on real data rather than only reading about.
	 */
	@Test
	public void segmentsOfOneVirusCountSeparately() throws IOException {
		assertEquals(8, admitAll("viral_segments.fna"));
	}

	/** So a selection below the segment count takes a virus apart, which a caller should know. */
	@Test
	public void aSelectionBelowTheSegmentCountTruncatesTheVirus() throws IOException {
		assertEquals(3, read("viral_segments.fna",
				selectionOf("NC_029549.1", "NC_029550.1", "NC_029551.1")).included.size());
	}

	// ---- harness --------------------------------------------------------------------------------

	/** Feeds every accession of the excerpt through the selection and returns how many genomes it took. */
	private int admitAll(String resource) throws IOException {
		GenomeKeyTrie selection = new GenomeKeyTrie();
		int genomes = 0;
		for (String line : Files.readAllLines(resourceFile(resource).toPath())) {
			if (line.startsWith(">")) {
				byte[] accession = line.substring(1, line.indexOf(' ')).getBytes(StandardCharsets.US_ASCII);
				if (selection.admit(accession, 0, accession.length)) {
					genomes++;
				}
			}
		}
		return genomes;
	}

	private static GenomeKeyTrie selectionOf(String... accessions) {
		GenomeKeyTrie selection = new GenomeKeyTrie();
		for (String accession : accessions) {
			byte[] bytes = accession.getBytes(StandardCharsets.US_ASCII);
			selection.admit(bytes, 0, bytes.length);
		}
		return selection;
	}

	private File resourceFile(String resource) {
		return new File(getClass().getResource("/refseq/" + resource).getFile());
	}

	private Reader read(String resource, GenomeKeyTrie selection) throws IOException {
		Reader reader = new Reader(tree, selection);
		reader.readFasta(resourceFile(resource));
		return reader;
	}

	private static void write(File file, String... lines) throws IOException {
		try (PrintWriter pw = new PrintWriter(file, StandardCharsets.UTF_8.name())) {
			for (String line : lines) {
				pw.println(line);
			}
		}
	}

	/** The tax ids the RefSeq catalog gives these accessions, keyed by their WGS or replicon prefix. */
	private static final class CatalogAccessionMap implements AccessionMap {
		private final TaxTree tree;
		private final GenomeKeyTrie admittedGenomes;

		private CatalogAccessionMap(TaxTree tree, GenomeKeyTrie admittedGenomes) {
			this.tree = tree;
			this.admittedGenomes = admittedGenomes;
		}

		@Override
		public void put(byte[] array, int start, int end, TaxIdNode node) {
		}

		@Override
		public TaxIdNode get(byte[] array, int start, int end, boolean assemblyAccessionsOnly) {
			String accession = new String(array, start, end - start, StandardCharsets.US_ASCII);
			if (accession.startsWith("NZ_JARSVH")) {
				return tree.getNodeByTaxId("1428");
			}
			if (accession.startsWith("NZ_JARSRC")) {
				return tree.getNodeByTaxId("1423");
			}
			if (accession.startsWith("NC_0295")) {
				return tree.getNodeByTaxId("1980428");
			}
			return null;
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

	/** Records which contigs were let in. */
	private static final class Reader extends AbstractRefSeqFastaReader {
		private final List<String> included = new ArrayList<>();

		private Reader(TaxTree tree, GenomeKeyTrie selection) {
			super(4096, allNodes(tree), new CatalogAccessionMap(tree, selection), 31, 1, false,
					new StringLong2DigitTrie());
		}

		private static Set<TaxIdNode> allNodes(TaxTree tree) {
			Set<TaxIdNode> nodes = new HashSet<>();
			for (String taxId : new String[] { "1423", "1428", "1980428" }) {
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
		}
	}
}
