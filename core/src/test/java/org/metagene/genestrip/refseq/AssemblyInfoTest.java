/*
 * Genestrip
 */
package org.metagene.genestrip.refseq;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertSame;
import static org.junit.Assert.assertTrue;

import org.junit.Test;
import org.metagene.genestrip.genbank.AssemblySummaryReader.AssemblyQuality;

/**
 * Tests the record the assembly metadata pass stores per genome.
 * <p>
 * The level is kept rather than consumed by the admit test, so that a build which filters nothing can
 * still be told which of its genomes are complete. That is what a topology inferred from complete
 * genomes alone needs, and it is the reason this class holds two things instead of one.
 */
public class AssemblyInfoTest {

	private static AssemblyInfo of(AssemblyQuality level) {
		return new AssemblyInfo(new byte[] { 'N', 'C', '_', '0', '0', '1' }, level);
	}

	/** The key is handed back as given: it is what every sequence of the assembly is filed under. */
	@Test
	public void testKeyIsKept() {
		byte[] key = new byte[] { 'N', 'C', '_', '9' };
		assertArrayEquals(key, new AssemblyInfo(key, AssemblyQuality.COMPLETE_LATEST).getKey());
	}

	/** Complete is complete whether or not it is the latest version of the assembly. */
	@Test
	public void testCompleteIsRecognised() {
		for (AssemblyQuality q : new AssemblyQuality[] { AssemblyQuality.COMPLETE_LATEST,
				AssemblyQuality.COMPLETE }) {
			assertTrue(q + " must count as complete", of(q).isComplete());
			assertFalse(q + " must not count as chromosome level", of(q).isChromosome());
		}
	}

	/** Chromosome level is not complete, which is the distinction a backbone turns on. */
	@Test
	public void testChromosomeIsNotComplete() {
		for (AssemblyQuality q : new AssemblyQuality[] { AssemblyQuality.CHROMOSOME_LATEST,
				AssemblyQuality.CHROMOSOME }) {
			assertTrue(q + " must count as chromosome level", of(q).isChromosome());
			assertFalse(q + " must not count as complete", of(q).isComplete());
		}
	}

	/** A draft is neither, so a consumer asking either question gets a plain no. */
	@Test
	public void testDraftIsNeither() {
		for (AssemblyQuality q : new AssemblyQuality[] { AssemblyQuality.SCAFFOLD_LATEST,
				AssemblyQuality.SCAFFOLD, AssemblyQuality.CONTIG_LATEST, AssemblyQuality.CONTIG }) {
			assertFalse(q + " must not count as complete", of(q).isComplete());
			assertFalse(q + " must not count as chromosome level", of(q).isChromosome());
		}
	}

	/** The level itself is available, for a consumer that wants more than the two questions. */
	@Test
	public void testLevelIsKept() {
		assertSame(AssemblyQuality.CHROMOSOME_LATEST, of(AssemblyQuality.CHROMOSOME_LATEST).getLevel());
	}
}
