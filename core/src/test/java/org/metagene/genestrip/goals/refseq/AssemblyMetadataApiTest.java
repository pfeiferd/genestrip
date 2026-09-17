/*
 * Genestrip
 */
package org.metagene.genestrip.goals.refseq;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.io.File;

import org.junit.Test;
import org.metagene.genestrip.APITest;
import org.metagene.genestrip.GSCommon;
import org.metagene.genestrip.GSGoalKey;
import org.metagene.genestrip.GSMaker;
import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.refseq.AccessionTrie;
import org.metagene.genestrip.refseq.AssemblyInfo;
import org.metagene.genestrip.refseq.GenomeKeyTrie;

/**
 * Tests that the assembly metadata is reachable the way a caller outside this package would reach it.
 * <p>
 * The value of the goal is a trie keyed by the genome key rather than by an accession, so a caller
 * holding {@code NC_001911.1} cannot look it up without deriving that key first. These tests pin the
 * accessors that do the deriving, since they are the whole of the interface a consumer of this goal
 * has to learn.
 */
public class AssemblyMetadataApiTest {

	@SuppressWarnings("unchecked")
	private AssemblyMetadataGoal<GSProject> goal() throws Exception {
		GSProject project = new GSProject(new GSCommon(APITest.getBaseDir()), "human_virus", true);
		return (AssemblyMetadataGoal<GSProject>) new GSMaker<GSProject>(project)
				.getGoal(GSGoalKey.ASSEMBLYMETA);
	}

	/** The goal exposes the two accessors a consumer needs, and they are public. */
	@Test
	public void testAccessorsArePublic() throws Exception {
		assertNotNull(AssemblyMetadataGoal.class.getMethod("getAssemblyInfo", String.class));
		assertNotNull(AssemblyMetadataGoal.class.getMethod("getAssemblyInfo", byte[].class, int.class,
				int.class));
		assertNotNull(goal());
	}

	/**
	 * The key an accession is filed under is not the accession. This is what the accessors hide, and
	 * the reason a caller cannot simply hand the trie a name.
	 */
	@Test
	public void testTheKeyIsNotTheAccession() {
		byte[] finished = key("NC_001911.1");
		assertEquals("the version is cut from a finished replicon", "NC_001911",
				new String(finished, 0, GenomeKeyTrie.genomeKeyLength(finished, 0, finished.length)));
		byte[] wgs = key("NZ_CABEIU010000001.1");
		assertEquals("a shotgun contig is cut back to its project prefix", "NZ_CABEIU",
				new String(wgs, 0, GenomeKeyTrie.genomeKeyLength(wgs, 0, wgs.length)));
	}

	/**
	 * Every sequence of one assembly resolves to the same record, which is what lets a caller count
	 * genomes rather than replicons.
	 */
	@Test
	public void testSequencesOfOneAssemblyShareARecord() {
		AccessionTrie<AssemblyInfo> trie = new AccessionTrie<AssemblyInfo>();
		AssemblyInfo info = new AssemblyInfo(key("NC_000913"),
				org.metagene.genestrip.genbank.AssemblySummaryReader.AssemblyQuality.COMPLETE_LATEST, false);
		for (String acc : new String[] { "NC_000913", "NC_000914" }) {
			byte[] b = key(acc);
			trie.set(b, 0, b.length, info);
		}
		for (String acc : new String[] { "NC_000913", "NC_000914" }) {
			byte[] b = key(acc);
			AssemblyInfo found = trie.get(b, 0, b.length);
			assertNotNull(acc + " must resolve", found);
			assertTrue(acc + " must be reported complete", found.isComplete());
			assertEquals("and to the same assembly", "NC_000913",
					new String(found.getKey()));
		}
		byte[] absent = key("NC_999999");
		assertNull("an accession of no known assembly must be null", trie.get(absent, 0, absent.length));
	}

	private static byte[] key(String s) {
		byte[] b = new byte[s.length()];
		for (int i = 0; i < s.length(); i++) {
			b[i] = (byte) s.charAt(i);
		}
		return b;
	}
}
