/*
 * Genestrip
 */
package org.metagene.genestrip.refseq;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.io.File;

import org.junit.BeforeClass;
import org.junit.Test;
import org.metagene.genestrip.genbank.AssemblySummaryReader.AssemblyQuality;

/**
 * Tests the assembly summary index against the real summary file where one is installed.
 * <p>
 * The file is large and is not in the repository, so every test here is skipped rather than failed
 * when it is absent: a build on a machine without a downloaded RefSeq must not go red for that. Set
 * {@code -Dgenestrip.refseqDir} to point at the directory holding it.
 */
public class AssemblySizeIndexTest {

	private static AssemblySizeIndex index;

	@BeforeClass
	public static void setUp() throws Exception {
		File f = summaryFile();
		if (f != null) {
			index = new AssemblySizeIndex(f);
		}
	}

	/** Returns the summary file, or null where none is installed. */
	private static File summaryFile() {
		String dir = System.getProperty("genestrip.refseqDir");
		if (dir == null) {
			dir = "../ft-db-exp2/data/common/refseq";
		}
		File f = new File(dir, "assembly_summary_refseq.txt");
		return f.exists() ? f : null;
	}

	/** Skips a test where no summary is installed. */
	private boolean skip() {
		return index == null;
	}

	@Test
	public void testIndexIsBuiltAndSorted() {
		if (skip()) {
			return;
		}
		// Every assembly is indexed under its taxon and its species and under both size columns, so
		// there are several entries per assembly and a real RefSeq has a lot of them.
		assertTrue("expected a substantial index, got " + index.size(), index.size() > 100000);
	}

	@Test
	public void testFindsAKnownCompleteGenome() {
		if (skip()) {
			return;
		}
		// M. tuberculosis: its chromosome is about 4.41 Mbp and the species has hundreds of closed
		// assemblies, so an exact hit must resolve to a complete one.
		AssemblySizeIndex.Match m = index.nearest(1773, 4411532, 0);
		assertNotNull("expected to locate the M. tuberculosis chromosome", m);
		assertEquals(0, m.getDistance());
		assertTrue("expected a complete assembly, got " + m.getLevel(),
				m.getLevel() == AssemblyQuality.COMPLETE_LATEST || m.getLevel() == AssemblyQuality.COMPLETE);
	}

	@Test
	public void testGapAdmitsANearMissAndExactDoesNot() {
		if (skip()) {
			return;
		}
		long size = 4411532;
		assertNotNull(index.nearest(1773, size, 0));
		// Shifted by more than the gap allows, an exact search must find nothing while a wide one does.
		assertNull("a size 50 off must not match at gap 0", index.nearest(1773, size + 50, 0));
		assertNotNull("a size 50 off must match at gap 100", index.nearest(1773, size + 50, 100));
	}

	@Test
	public void testNearestWinsWhereSeveralAreInRange() {
		if (skip()) {
			return;
		}
		// With a wide gap several assemblies of a well sequenced species fall in range; whichever is
		// returned must be the nearest, so its distance cannot exceed that of an exact neighbour.
		AssemblySizeIndex.Match wide = index.nearest(1280, 2821361, 10000);
		assertNotNull(wide);
		assertTrue("distance must be within the gap", wide.getDistance() <= 10000);
		AssemblySizeIndex.Match narrow = index.nearest(1280, 2821361, 10);
		if (narrow != null) {
			assertTrue("a wider search must not return a worse match",
					wide.getDistance() <= narrow.getDistance());
		}
	}

	@Test
	public void testUnknownTaxonAndImpossibleSizeFindNothing() {
		if (skip()) {
			return;
		}
		assertNull("no assembly of a made-up taxon", index.nearest(999999999, 4411532, 100));
		assertNull("no assembly of one base", index.nearest(1773, 1, 100));
	}

	/**
	 * The gap is a tolerance in both directions. Only the positive side was covered before, which
	 * would not have caught a comparison written one-sided - and the sum of a run of replicons falls
	 * short of the recorded length as readily as it overshoots it.
	 * <p>
	 * The assertion deliberately does not claim that a near miss finds nothing: a taxon as densely
	 * sequenced as this one has assemblies a couple of bases apart, so any such claim would be about
	 * RefSeq's contents rather than about the lookup.
	 */
	@Test
	public void testGapAppliesBelowAsWellAsAbove() {
		if (skip()) {
			return;
		}
		long size = 4411532;
		assertNotNull("exact match expected", index.nearest(1773, size, 0));
		for (int d : new int[] { 1, 50, 100, 1000 }) {
			assertNotNull(d + " below must match at gap " + d, index.nearest(1773, size - d, d));
			assertNotNull(d + " above must match at gap " + d, index.nearest(1773, size + d, d));
		}
	}

	/** A match reports how far off it was as an unsigned distance, never wider than the gap. */
	@Test
	public void testMatchDistanceIsUnsignedAndWithinGap() {
		if (skip()) {
			return;
		}
		long size = 4411532;
		assertEquals("an exact hit is at distance 0", 0, index.nearest(1773, size, 0).getDistance());
		for (int d : new int[] { 1, 50, 100, 1000 }) {
			long below = index.nearest(1773, size - d, d).getDistance();
			long above = index.nearest(1773, size + d, d).getDistance();
			assertTrue("distance below must not be negative, was " + below, below >= 0);
			assertTrue("distance above must not be negative, was " + above, above >= 0);
			assertTrue("distance below must be within the gap, was " + below, below <= d);
			assertTrue("distance above must be within the gap, was " + above, above <= d);
		}
	}

	/**
	 * Widening the gap must never return a worse match. The scan keeps the minimum over the whole
	 * window rather than stopping at the first entry in range, so the distance is monotonically
	 * non-increasing in the gap and settles once the true nearest assembly is inside it.
	 */
	@Test
	public void testBestMatchHoldsAsTheGapWidens() {
		if (skip()) {
			return;
		}
		int[] taxids = { 1773, 1280, 562, 1313 };
		long[] sizes = { 4411532, 2821361, 4641652, 2038615 };
		for (int t = 0; t < taxids.length; t++) {
			long previous = Long.MAX_VALUE;
			for (int gap : new int[] { 0, 10, 100, 1000, 10000, 1000000, 100000000 }) {
				AssemblySizeIndex.Match m = index.nearest(taxids[t], sizes[t] + 7777, gap);
				if (m == null) {
					continue;
				}
				assertTrue("distance must stay within the gap for taxon " + taxids[t],
						m.getDistance() <= gap);
				assertTrue("widening the gap worsened the match for taxon " + taxids[t] + ": "
						+ previous + " -> " + m.getDistance(), m.getDistance() <= previous);
				previous = m.getDistance();
			}
			assertTrue("a wide search must find something for taxon " + taxids[t],
					previous < Long.MAX_VALUE);
		}
	}

	/**
	 * A gap wide enough to cover every assembly of a taxon must return the same match as a gap wide
	 * enough to cover the whole size range, which is what shows the scan is not truncating.
	 */
	@Test
	public void testAVeryWideGapIsStable() {
		if (skip()) {
			return;
		}
		AssemblySizeIndex.Match wide = index.nearest(1280, 2821361 + 7777, 100000000);
		AssemblySizeIndex.Match wider = index.nearest(1280, 2821361 + 7777, Integer.MAX_VALUE);
		assertNotNull(wide);
		assertNotNull(wider);
		assertEquals("the match must not change once every assembly is in range",
				wide.getDistance(), wider.getDistance());
		assertEquals(wide.getLevel(), wider.getLevel());
	}

	/**
	 * The index holds admitted levels only, so whatever a search returns is of an admitted level.
	 * That is what lets the caller take a hit as an answer about completeness instead of screening
	 * the level afterwards.
	 */
	@Test
	public void testOnlyAdmittedLevelsAreIndexed() {
		if (skip()) {
			return;
		}
		int[] taxids = { 1773, 1280, 562, 1313, 287, 573 };
		int found = 0;
		for (int taxid : taxids) {
			AssemblySizeIndex.Match m = index.nearest(taxid, 3000000, Integer.MAX_VALUE);
			if (m == null) {
				continue;
			}
			found++;
			assertTrue("a match must be of an admitted level, was " + m.getLevel(),
					AssemblySizeIndex.COMPLETE_OR_CHROMOSOME.contains(m.getLevel()));
		}
		assertTrue("expected a match for at least one well sequenced taxon", found > 0);
	}

	/**
	 * A draft nearer in size than the admitted assembly must not shadow it. Reading the same summary
	 * over all levels is what the old behaviour amounted to, and it finds a contig assembly three
	 * bases off where the admitting index finds the assembly the caller is asking about.
	 */
	@Test
	public void testADraftDoesNotShadowAnAdmittedAssembly() throws Exception {
		if (skip()) {
			return;
		}
		AssemblySizeIndex all = new AssemblySizeIndex(summaryFile(), AssemblySizeIndex.ALL_LEVELS);
		AssemblySizeIndex.Match loose = all.nearest(1280, 2829138, 100);
		AssemblySizeIndex.Match strict = index.nearest(1280, 2829138, 100);
		assertNotNull("the all-levels index finds a draft here", loose);
		assertNotNull("the admitting index must still find the assembly it is asked about", strict);
		assertTrue("the all-levels match is the draft that used to win",
				loose.getDistance() < strict.getDistance());
		assertTrue("the all-levels match is not of an admitted level",
				!AssemblySizeIndex.COMPLETE_OR_CHROMOSOME.contains(loose.getLevel()));
		assertTrue("and the admitting index's match is",
				AssemblySizeIndex.COMPLETE_OR_CHROMOSOME.contains(strict.getLevel()));
	}

	/** The all-levels index is available for a caller wanting the summary's metadata, and is larger. */
	@Test
	public void testAllLevelsIndexesMore() throws Exception {
		if (skip()) {
			return;
		}
		AssemblySizeIndex all = new AssemblySizeIndex(summaryFile(), AssemblySizeIndex.ALL_LEVELS);
		assertTrue("indexing every level must yield more entries: " + all.size() + " vs " + index.size(),
				all.size() > index.size() * 4);
	}

	/**
	 * The two-neighbour search must return what a scan of the whole window would have returned. The
	 * entries of a taxon are sorted by size, so the nearest to a sum is one of its two neighbours and
	 * the window between them holds nothing nearer - this checks that over many sums and gaps rather
	 * than trusting the argument.
	 */
	@Test
	public void testTwoNeighbourSearchAgreesWithAFullScan() {
		if (skip()) {
			return;
		}
		int[] taxids = { 1773, 1280, 562, 1313, 287, 573, 1496, 90370 };
		int compared = 0;
		for (int taxid : taxids) {
			for (long base : new long[] { 1000000, 2829138, 3000000, 4411532, 5000000 }) {
				for (int off : new int[] { -5000, -77, 0, 77, 5000 }) {
					for (int gap : new int[] { 0, 1, 100, 10000, 1000000, Integer.MAX_VALUE }) {
						long sum = base + off;
						AssemblySizeIndex.Match m = index.nearest(taxid, sum, gap);
						Long scanned = index.nearestByScan(taxid, sum, gap);
						if (scanned == null) {
							assertNull("scan found nothing but the search did, taxon " + taxid + " sum " + sum
									+ " gap " + gap, m);
						} else {
							assertNotNull("scan found a match but the search did not, taxon " + taxid + " sum "
									+ sum + " gap " + gap, m);
							assertEquals("distance differs, taxon " + taxid + " sum " + sum + " gap " + gap,
									(long) scanned, m.getDistance());
						}
						compared++;
					}
				}
			}
		}
		assertEquals("every combination must have been compared", 8 * 5 * 5 * 6, compared);
	}

	/**
	 * Building over complete assemblies alone must exclude chromosome-level ones, which is what the
	 * {@code complete} setting asks for. The levels are decided when the index is built rather than
	 * screened after a search, so this is where the setting takes effect.
	 */
	@Test
	public void testCompleteOnlyExcludesChromosomeLevel() throws Exception {
		if (skip()) {
			return;
		}
		AssemblySizeIndex completeOnly = new AssemblySizeIndex(summaryFile(), AssemblySizeIndex.COMPLETE_ONLY);
		assertTrue("admitting chromosome level must index more: " + index.size() + " vs " + completeOnly.size(),
				index.size() > completeOnly.size());
		int checked = 0;
		for (int taxid : new int[] { 1773, 1280, 562, 1313, 287, 573 }) {
			AssemblySizeIndex.Match m = completeOnly.nearest(taxid, 3000000, Integer.MAX_VALUE);
			if (m == null) {
				continue;
			}
			checked++;
			assertTrue("a match must be complete, was " + m.getLevel(), m.isComplete());
			assertFalse("and must not be chromosome level", m.isChromosome());
		}
		assertTrue("expected a match for at least one taxon", checked > 0);
	}

	/**
	 * The two flags say which of the admitted levels was found, and say it exclusively: a caller
	 * reckoning with how much of an organism the database holds needs to tell a gapless assembly from
	 * one carrying about a per cent of gaps.
	 */
	@Test
	public void testCompleteAndChromosomeFlagsAreExclusive() {
		if (skip()) {
			return;
		}
		int complete = 0;
		int chromosome = 0;
		for (int taxid = 2; taxid < 60000; taxid++) {
			AssemblySizeIndex.Match m = index.nearest(taxid, 3000000, Integer.MAX_VALUE);
			if (m == null) {
				continue;
			}
			assertTrue("a match must be one of the admitted levels, was " + m.getLevel(),
					m.isComplete() != m.isChromosome());
			if (m.isComplete()) {
				complete++;
			} else {
				chromosome++;
			}
		}
		assertTrue("expected complete matches", complete > 0);
		assertTrue("expected chromosome-level matches, else the flag is untested", chromosome > 0);
	}

	/**
	 * Every entry of one assembly carries one id, and different assemblies carry different ones. The
	 * id stands in for the assembly accession the accession catalog does not have, so a caller can
	 * tell a genome it has already accounted for from a new one.
	 */
	@Test
	public void testAssemblyIdsIdentifyAssemblies() {
		if (skip()) {
			return;
		}
		assertTrue("every entry must belong to an assembly",
				index.getAssemblyCount() > 0 && index.getAssemblyCount() <= index.size());
		// An assembly is indexed under its taxon and its species and under both size columns, so there
		// are several entries per assembly and never more assemblies than entries.
		assertTrue("entries must outnumber assemblies", index.size() > index.getAssemblyCount());
		int distinct = 0;
		java.util.Set<Integer> ids = new java.util.HashSet<Integer>();
		for (int taxid : new int[] { 1773, 1280, 562, 1313, 287, 573 }) {
			AssemblySizeIndex.Match m = index.nearest(taxid, 3000000, Integer.MAX_VALUE);
			if (m == null) {
				continue;
			}
			assertTrue("an id must be in range",
					m.getAssemblyId() >= 0 && m.getAssemblyId() < index.getAssemblyCount());
			if (ids.add(m.getAssemblyId())) {
				distinct++;
			}
		}
		assertTrue("different taxa must find different assemblies", distinct > 1);
	}

	/** The same query must name the same assembly every time, or a caller could not dedupe on it. */
	@Test
	public void testAssemblyIdIsStable() {
		if (skip()) {
			return;
		}
		AssemblySizeIndex.Match a = index.nearest(1280, 2829138, 100);
		AssemblySizeIndex.Match b = index.nearest(1280, 2829138, 100);
		assertNotNull(a);
		assertNotNull(b);
		assertEquals(a.getAssemblyId(), b.getAssemblyId());
	}
}
