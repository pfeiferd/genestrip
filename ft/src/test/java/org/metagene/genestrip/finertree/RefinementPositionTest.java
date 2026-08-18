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
package org.metagene.genestrip.finertree;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertSame;
import static org.junit.Assert.assertTrue;

import java.util.Arrays;
import java.util.Collections;

import org.junit.Test;
import org.metagene.genestrip.finertree.FTConfigKey.RefinementPosition;
import org.metagene.genestrip.finertree.FTConfigKey.RefinementPosition.Limit;
import org.metagene.genestrip.tax.Rank;
import org.metagene.genestrip.tax.TaxIdInfo;

/**
 * Tests which taxonomy nodes a configured refinement position selects.
 * <p>
 * This decides where a refinement happens at all: {@code kmerindexbloom} turns the answer into the
 * bit set its per-k-mer hot path tests, and a node wrongly selected has its k-mers pushed down to
 * nodes that should never have been created. The ranks are ordered from the root downwards, so a
 * rank <em>larger</em> than another is the one closer to the root and carries the smaller ordinal -
 * which is the point most easily got backwards, and hence the one tested here explicitly.
 */
public class RefinementPositionTest {
	/** A node of a given rank and tax id, with an optional parent. */
	private static class Node extends TaxIdInfo {
		private static final long serialVersionUID = 1L;

		private final TaxIdInfo parent;

		Node(String taxId, Rank rank, TaxIdInfo parent) {
			super(taxId, rank);
			this.parent = parent;
		}

		@Override
		public TaxIdInfo getParent() {
			return parent;
		}
	}

	private static Node node(String taxId, Rank rank) {
		return new Node(taxId, rank, null);
	}

	/**
	 * The comparison of a configured position is read off the front of its token, so that the rank
	 * may follow it directly. Two-character operators must win over the one-character ones they begin
	 * with, or {@code >=genus} would be read as {@code >} and quietly select one rank too few.
	 */
	@Test
	public void testLimitsAreParsedFromTheFrontOfTheirTokens() {
		assertSame(Limit.EQ, Limit.fromString("="));
		assertSame(Limit.LARGER, Limit.fromString(">"));
		assertSame(Limit.LESS, Limit.fromString("<"));
		assertSame(Limit.LARGER_EQ, Limit.fromString(">="));
		assertSame(Limit.LESS_EQ, Limit.fromString("<="));
		assertSame(Limit.ALL, Limit.fromString("*"));

		assertSame("the longer operator must win over its own prefix", Limit.LARGER_EQ,
				Limit.fromString(">=genus"));
		assertSame(Limit.LESS_EQ, Limit.fromString("<=species"));
		assertSame(Limit.LARGER, Limit.fromString(">genus"));
		assertSame(Limit.LESS, Limit.fromString("<species"));
		assertSame(Limit.EQ, Limit.fromString("=genus"));

		assertNull("a token starting with none of the operators has no limit", Limit.fromString("genus"));
		assertNull(Limit.fromString(""));
		for (Limit limit : Limit.values()) {
			assertSame("every limit should parse back from its own token", limit,
					Limit.fromString(limit.getComp()));
		}
	}

	/** The sentinel of the default configuration selects every node there is. */
	@Test
	public void testAllPositionsMatchesEverything() {
		for (Rank rank : new Rank[] { Rank.SUPERKINGDOM, Rank.GENUS, Rank.SPECIES, Rank.STRAIN }) {
			assertTrue(rank + " should be matched by the catch-all position",
					RefinementPosition.ALL_POSITIONS.isMatchingNodeForPosition(node("1", rank)));
		}
	}

	/** A position naming a rank selects that rank and no other. */
	@Test
	public void testEqualRankMatchesThatRankAlone() {
		RefinementPosition genus = new RefinementPosition(Rank.GENUS, null, Limit.EQ);
		assertTrue(genus.isMatchingNodeForPosition(node("1", Rank.GENUS)));
		assertFalse(genus.isMatchingNodeForPosition(node("1", Rank.SPECIES)));
		assertFalse(genus.isMatchingNodeForPosition(node("1", Rank.FAMILY)));
	}

	/** A position naming a tax id selects that node and no other, whatever its rank. */
	@Test
	public void testEqualTaxIdMatchesThatNodeAlone() {
		RefinementPosition byTaxId = new RefinementPosition(null, "1496", Limit.EQ);
		assertTrue(byTaxId.isMatchingNodeForPosition(node("1496", Rank.SPECIES)));
		assertFalse(byTaxId.isMatchingNodeForPosition(node("1497", Rank.SPECIES)));
	}

	/**
	 * "Larger than a rank" means closer to the root. Ranks are ordered downwards, so this is the
	 * comparison that reads the wrong way round if one goes by the ordinals alone.
	 */
	@Test
	public void testLargerSelectsTheRanksAboveIt() {
		RefinementPosition aboveGenus = new RefinementPosition(Rank.GENUS, null, Limit.LARGER);
		assertTrue("family is above genus", aboveGenus.isMatchingNodeForPosition(node("1", Rank.FAMILY)));
		assertTrue("phylum is above genus", aboveGenus.isMatchingNodeForPosition(node("1", Rank.PHYLUM)));
		assertFalse("genus itself is not above genus", aboveGenus.isMatchingNodeForPosition(node("1", Rank.GENUS)));
		assertFalse("species is below genus", aboveGenus.isMatchingNodeForPosition(node("1", Rank.SPECIES)));

		RefinementPosition genusAndAbove = new RefinementPosition(Rank.GENUS, null, Limit.LARGER_EQ);
		assertTrue("the rank itself is included", genusAndAbove.isMatchingNodeForPosition(node("1", Rank.GENUS)));
		assertTrue(genusAndAbove.isMatchingNodeForPosition(node("1", Rank.FAMILY)));
		assertFalse(genusAndAbove.isMatchingNodeForPosition(node("1", Rank.SPECIES)));
	}

	/** "Less than a rank" means further from the root, i.e. the more specific ranks. */
	@Test
	public void testLessSelectsTheRanksBelowIt() {
		RefinementPosition belowGenus = new RefinementPosition(Rank.GENUS, null, Limit.LESS);
		assertTrue("species is below genus", belowGenus.isMatchingNodeForPosition(node("1", Rank.SPECIES)));
		assertFalse("genus itself is not below genus", belowGenus.isMatchingNodeForPosition(node("1", Rank.GENUS)));
		assertFalse("family is above genus", belowGenus.isMatchingNodeForPosition(node("1", Rank.FAMILY)));

		RefinementPosition genusAndBelow = new RefinementPosition(Rank.GENUS, null, Limit.LESS_EQ);
		assertTrue("the rank itself is included", genusAndBelow.isMatchingNodeForPosition(node("1", Rank.GENUS)));
		assertTrue(genusAndBelow.isMatchingNodeForPosition(node("1", Rank.SPECIES)));
		assertFalse(genusAndBelow.isMatchingNodeForPosition(node("1", Rank.FAMILY)));
	}

	/**
	 * A node of no determinate rank - which the taxonomy is full of - has no ordinal to compare, so
	 * the decision is taken by walking up to a node that has one.
	 */
	@Test
	public void testIndeterminateRankIsDecidedByItsAncestors() {
		Node family = new Node("100", Rank.FAMILY, null);
		Node genus = new Node("200", Rank.GENUS, family);
		Node noRankBelowGenus = new Node("300", Rank.NO_RANK, genus);
		Node noRankBelowFamily = new Node("400", Rank.NO_RANK, family);

		RefinementPosition belowGenus = new RefinementPosition(Rank.GENUS, null, Limit.LESS);
		assertTrue("a rankless node under a genus counts as below it",
				belowGenus.isMatchingNodeForPosition(noRankBelowGenus));
		assertFalse("a rankless node that has no genus above it does not",
				belowGenus.isMatchingNodeForPosition(noRankBelowFamily));

		RefinementPosition aboveGenus = new RefinementPosition(Rank.GENUS, null, Limit.LARGER);
		assertFalse(aboveGenus.isMatchingNodeForPosition(noRankBelowGenus));
		assertTrue(aboveGenus.isMatchingNodeForPosition(noRankBelowFamily));
	}

	/** Of several configured positions the first matching one is returned, and none means null. */
	@Test
	public void testFirstMatchingPositionWins() {
		RefinementPosition species = new RefinementPosition(Rank.SPECIES, null, Limit.EQ);
		RefinementPosition genus = new RefinementPosition(Rank.GENUS, null, Limit.EQ);
		assertSame(genus, RefinementPosition.getMatchingNodeFor(node("1", Rank.GENUS),
				Arrays.asList(species, genus)));
		assertSame(species, RefinementPosition.getMatchingNodeFor(node("1", Rank.SPECIES),
				Arrays.asList(species, genus)));
		assertNull("a family matches neither", RefinementPosition.getMatchingNodeFor(node("1", Rank.FAMILY),
				Arrays.asList(species, genus)));
		assertNull("nothing configured matches nothing", RefinementPosition.getMatchingNodeFor(
				node("1", Rank.GENUS), Collections.emptyList()));
	}

	/** The catch-all in a list makes every node match, whatever stands beside it. */
	@Test
	public void testCatchAllInAListMatchesAnyNode() {
		assertEquals(RefinementPosition.ALL_POSITIONS, RefinementPosition.getMatchingNodeFor(
				node("1", Rank.SUBSPECIES),
				Arrays.asList(RefinementPosition.ALL_POSITIONS,
						new RefinementPosition(Rank.GENUS, null, Limit.EQ))));
	}
}
