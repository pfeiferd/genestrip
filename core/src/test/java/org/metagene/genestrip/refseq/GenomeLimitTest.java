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
import static org.junit.Assert.assertNotNull;
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

	/**
	 * What a limit over genomes counts. Only the accession kinds an assembly is made of do, and that
	 * is independent of {@code refseq.assemblyAccessionsOnly}: a region record is a few kilobases out
	 * of a genome and never an assembly of one, so it must take no genome's place whether or not its
	 * sequence enters the database. Counting it did cost the pneumococcus of the paper's `strepto'
	 * twelve of its fifty places, against 9,263 assemblies on offer -- and the genome key cannot
	 * catch that, since a region record stands for itself exactly as a chromosome does.
	 */
	@Test
	public void onlyAssemblyAccessionsCountAsGenomes() {
		for (String assembly : new String[] { "NC_003028.3", "NZ_CABEIU010000001.1", "NZ_AP017971.1",
				"AC_000091.1" }) {
			assertTrue(assembly + " is what an assembly is made of", isAssembly(assembly));
		}
		for (String region : new String[] { "NG_048006.1", "NT_187512.1", "NW_025791637.1" }) {
			assertFalse(region + " is a region of a genome, not one", isAssembly(region));
		}
	}

	/** And the key does not help here: both kinds stand for themselves, so both get their own. */
	@Test
	public void aRegionRecordIsNotToldFromAChromosomeByItsKey() {
		assertEquals("NG_048006", key("NG_048006.1"));
		assertEquals("NC_003028", key("NC_003028.3"));
	}

	// ---- which genome a taxon is represented by ---------------------------------------------------

	/**
	 * The summary states a marked assembly under its own taxon and under its species, and both are
	 * worth a place: which of the two a limit is counted at depends on {@code maxPerTaxidRank}.
	 */
	@Test
	public void aMarkedAssemblyIsKnownByTaxonAndBySpecies() throws IOException {
		ReferenceGenomes refs = referenceGenomes(
				summaryRow("GCF_000000001.1", "na", "reference genome", "170187", "1313"),
				summaryRow("GCF_000000002.1", "na", "na", "1309", "1309"));
		assertTrue(refs.hasReference("170187"));
		assertTrue(refs.hasReference("1313"));
		assertFalse(refs.hasReference("1309"));
		assertEquals(2, refs.taxonCount());
	}

	/**
	 * A marked draft is recognised by the prefix its contigs carry. The summary states the WGS master
	 * in its GenBank form and the catalog carries the RefSeq one, so the key is cut from {@code NZ_}
	 * plus that form -- checked against one release, this finds 88 of the 89 Streptococcus species
	 * whose reference is a draft.
	 */
	@Test
	public void aMarkedDraftIsRecognisedByItsContigs() throws IOException {
		ReferenceGenomes refs = referenceGenomes(
				summaryRow("GCF_000000003.1", "AEVD00000000.1", "reference genome", "889204", "68892"));
		assertEquals(1, refs.draftCount());
		assertTrue(admitted(refs, "NZ_AEVD01000001.1"));
		assertTrue(admitted(refs, "NZ_AEVD01000273.1"));
		assertFalse(admitted(refs, "NZ_AEVF01000001.1"));
	}

	/** A finished assembly has no WGS master and is left to the length index, not to the prefix. */
	@Test
	public void aFinishedReferenceIsNotRecognisedByAPrefix() throws IOException {
		ReferenceGenomes refs = referenceGenomes(
				summaryRow("GCF_000000004.1", "na", "reference genome", "1313", "1313"));
		assertEquals(0, refs.draftCount());
		assertTrue(refs.hasReference("1313"));
	}

	/**
	 * The place kept for a marked assembly. It is held back only while one is expected and has not
	 * arrived, so a taxon without a marked assembly fills to the limit as it always did.
	 */
	@Test
	public void aPlaceIsKeptForTheMarkedAssembly() {
		// A taxon that has one, before it arrives: the last place is not given away.
		assertEquals(49, ReferenceGenomes.roomFor(50, false, false, true));
		// The marked assembly itself may take it.
		assertEquals(50, ReferenceGenomes.roomFor(50, true, false, true));
		// Once it is in, nothing is held back any more.
		assertEquals(50, ReferenceGenomes.roomFor(50, false, true, true));
		// And a taxon the summary marks nothing for is unaffected.
		assertEquals(50, ReferenceGenomes.roomFor(50, false, false, false));
		// Without a limit there is nothing to keep a place in.
		assertEquals(Integer.MAX_VALUE,
				ReferenceGenomes.roomFor(Integer.MAX_VALUE, false, false, true));
	}

	/**
	 * The rule over a sequence of genomes, as the catalog states them: the marked assembly gets in
	 * however late it arrives, and the limit is never exceeded. This is the loop
	 * {@code AccessionMapGoal.admit()} runs, with the two lookups of a real build -- which genome is
	 * marked, and how many the taxon has taken -- standing in for its trie and its counter.
	 */
	private static int[] admitInOrder(int maxGenomes, int genomes, int markedAt, boolean hasReference) {
		int taken = 0;
		boolean referenceAdmitted = false;
		int markedPosition = -1;
		for (int i = 0; i < genomes; i++) {
			boolean marked = hasReference && i == markedAt;
			if (taken < ReferenceGenomes.roomFor(maxGenomes, marked, referenceAdmitted, hasReference)) {
				taken++;
				if (marked) {
					referenceAdmitted = true;
					markedPosition = i;
				}
			}
		}
		return new int[] { taken, markedPosition };
	}

	/** The marked assembly is taken in even where the catalog states it after the limit is reached. */
	@Test
	public void theMarkedAssemblyGetsInHoweverLateItArrives() {
		int[] r = admitInOrder(50, 500, 499, true);
		assertEquals("the marked assembly must be in", 499, r[1]);
		assertEquals("and the limit must hold", 50, r[0]);
	}

	/** Holding a place back does not raise the limit: fifty stay fifty, one of them the marked one. */
	@Test
	public void keepingAPlaceDoesNotRaiseTheLimit() {
		assertEquals(50, admitInOrder(50, 500, 0, true)[0]);
		assertEquals(50, admitInOrder(50, 500, 49, true)[0]);
		assertEquals(50, admitInOrder(50, 500, 250, true)[0]);
	}

	/** A taxon the summary marks nothing for fills to the limit as it always did. */
	@Test
	public void ataxonWithoutAMarkedAssemblyIsUnaffected() {
		int[] r = admitInOrder(50, 500, -1, false);
		assertEquals(50, r[0]);
		assertEquals(-1, r[1]);
	}

	/**
	 * A taxon whose marked assembly the release does not carry ends one genome short. That is the
	 * price of holding the limit exactly instead of exceeding it by one, and it is worth pinning:
	 * silently handing the place to another genome would make the count depend on the catalog again.
	 */
	@Test
	public void aMarkedAssemblyTheReleaseLacksCostsOnePlace() {
		assertEquals(49, admitInOrder(50, 500, -1, true)[0]);
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
	 * A limit over genomes gates genomes. A region record and a transcript are not ones, are in no
	 * selection, and pass whatever the selection holds -- whether they enter at all is what
	 * {@code refseq.assemblyAccessionsOnly} and the sequence type decide. Letting the selection
	 * answer for them made a limit empty a database of its transcripts.
	 */
	@Test
	public void whatTheLimitDoesNotCountItDoesNotGate() throws IOException {
		List<String> included = readWithRegionAndRna(selectionOf("NZ_AAAA01000001.1")).included;
		assertTrue("the region record passes", included.contains("NG_048006.1"));
		assertTrue("and so does the transcript", included.contains("NR_000001.1"));
		assertFalse("while an assembly outside the selection stays out",
				included.contains("NZ_CP000001.1"));
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
	 * A genome outside the selection must not leave its taxon behind. The reader marks a contig's
	 * taxon as required and gives it its artificial data node, and {@code markRequired()} is the only
	 * thing that keeps a node in the {@link org.metagene.genestrip.tax.SmallTaxTree}. Asking about the
	 * selection afterwards - as this once did - left a node nothing is ever filed at standing in the
	 * database tree, where every measure that divides by the data taxa under a node counts it as a
	 * candidate that can never carry a k-mer.
	 */
	@Test
	public void aRejectedGenomeLeavesNoNodeBehind() throws IOException {
		readCreatingNodes(selectionOf("NZ_AAAA01000001.1"));

		TaxIdNode strain = tree.getNodeByTaxId("170187");
		assertFalse(strain.isRequired());
		assertNull(strain.getDataChild());
		assertNull(tree.toSmallTaxTree().getNodeByTaxId("170187"));

		// The admitted genome's own taxon is there with its data node, so the test cannot pass by
		// creating nothing at all.
		TaxIdNode species = tree.getNodeByTaxId("1313");
		assertTrue(species.isRequired());
		assertNotNull(species.getDataChild());
	}

	/** And a genome inside it gets its taxon and its data node, wherever the taxonomy files it. */
	@Test
	public void anAdmittedGenomeGetsItsNodes() throws IOException {
		readCreatingNodes(selectionOf("NZ_AAAA01000001.1", "NZ_BBBB01000001.1"));

		TaxIdNode strain = tree.getNodeByTaxId("170187");
		assertTrue(strain.isRequired());
		assertNotNull(strain.getDataChild());
		assertNotNull(tree.toSmallTaxTree().getNodeByTaxId("170187"));
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

	private static boolean isAssembly(String accession) {
		byte[] bytes = accession.getBytes(StandardCharsets.US_ASCII);
		return AccessionFileProcessor.isAssemblyAccession(bytes, 0);
	}

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

	/**
	 * Reads the same fasta as {@link #read(GenomeKeyTrie)}, but with a reader that creates the
	 * artificial nodes - which is what the sizing pass does, and the only pass that does.
	 */
	private void readCreatingNodes(GenomeKeyTrie selection) throws IOException {
		File fasta = folder.newFile();
		write(fasta,
				">NZ_AAAA01000001.1 first project, first contig", "ACGTACGTACGT",
				">NZ_BBBB01000001.1 second project, filed at the strain", "ACGTACGTACGT");
		new NodeCreatingReader(tree, selection).readFasta(fasta);
	}

	/** One row of an assembly summary, with the five columns {@link ReferenceGenomes} reads. */
	private static String summaryRow(String accession, String wgsMaster, String category, String taxId,
			String speciesTaxId) {
		StringBuilder row = new StringBuilder();
		row.append(accession).append('\t').append("PRJNA1").append('\t').append("SAMN1").append('\t')
				.append(wgsMaster).append('\t').append(category).append('\t').append(taxId).append('\t')
				.append(speciesTaxId).append('\t').append("an organism").append('\t').append("na");
		return row.toString();
	}

	private ReferenceGenomes referenceGenomes(String... rows) throws IOException {
		File summary = folder.newFile();
		String[] lines = new String[rows.length + 1];
		lines[0] = "#assembly_accession\tbioproject\tbiosample\twgs_master\trefseq_category\ttaxid"
				+ "\tspecies_taxid\torganism_name\tinfraspecific_name";
		System.arraycopy(rows, 0, lines, 1, rows.length);
		write(summary, lines);
		return new ReferenceGenomes(summary);
	}

	private static boolean admitted(ReferenceGenomes refs, String accession) {
		byte[] bytes = accession.getBytes(StandardCharsets.US_ASCII);
		return refs.isReferenceDraft(bytes, 0, bytes.length);
	}

	/**
	 * Reads a fasta of an admitted genome, a region record and a transcript against the given
	 * selection, which is what a limit has to leave alone.
	 */
	private Reader readWithRegionAndRna(GenomeKeyTrie selection) throws IOException {
		File fasta = folder.newFile();
		write(fasta,
				">NZ_AAAA01000001.1 first project, first contig", "ACGTACGTACGT",
				">NZ_CP000001.1 a finished replicon of another species", "ACGTACGTACGT",
				">NG_048006.1 a region of a genome, not an assembly of one", "ACGTACGTACGT",
				">NR_000001.1 a transcript", "ACGTACGTACGT");
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

	/** Creates the artificial data nodes for the contigs it lets in, as the sizing pass does. */
	private static final class NodeCreatingReader extends ReworkingStoreFastaReader {
		private NodeCreatingReader(TaxTree tree, GenomeKeyTrie selection) {
			super(4096, allNodes(tree), new LetterAccessionMap(tree, selection), 31, -1, 1, false,
					new StringLong2DigitTrie(), true, tree, true, false, false, false, true,
					counter -> "99" + counter, null);
		}

		private static Set<TaxIdNode> allNodes(TaxTree tree) {
			Set<TaxIdNode> nodes = new HashSet<>();
			for (String taxId : new String[] { "1313", "170187", "1309" }) {
				nodes.add(tree.getNodeByTaxId(taxId));
			}
			return nodes;
		}

		@Override
		protected boolean handleStore(long kmer) {
			return false;
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
