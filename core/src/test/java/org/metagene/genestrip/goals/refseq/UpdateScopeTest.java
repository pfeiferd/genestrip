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
package org.metagene.genestrip.goals.refseq;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.util.Collections;
import java.util.HashSet;
import java.util.Set;

import org.junit.Test;
import org.metagene.genestrip.GSConfigKey;
import org.metagene.genestrip.GSConfigKey.UpdateScope;
import org.metagene.genestrip.make.ConfigParamInfo;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;

/**
 * Tests {@link GSConfigKey#UPDATE_SCOPE} and the region selection it drives in
 * {@link DBGoal#isRegionInScope}.
 * <p>
 * The three scopes differ only in which regions of the RefSeq release the update is shown, so that
 * decision is what these tests pin down. The case that motivates {@link UpdateScope#OTHER_TAXA_ONLY}
 * is a database whose leaves are the genomes of the requested taxon themselves: the release then
 * holds those very genomes a second time under the taxon's own node, and letting the update see them
 * would raise every k-mer to that node. The project's own fastas are a different matter and are
 * covered by {@link #projectFastasAreNeverOutOfScope()}.
 */
public class UpdateScopeTest {
	/** A region of the RefSeq release, i.e. one the scope may restrict. */
	private static final boolean FROM_RELEASE = true;
	/** A region of a fasta the project supplies itself - {@code additional.txt} or Genbank. */
	private static final boolean FROM_PROJECT_FASTA = false;

	/** A requested taxon, a strain below it (so also in the selection), and one of another taxon. */
	private final TaxIdNode requested = new TaxIdNode("1496");
	private final TaxIdNode strainBelowRequested = new TaxIdNode("1163671");
	private final TaxIdNode other = new TaxIdNode("1280");
	/** A taxon `taxids.txt' struck out with a leading `-', and one below it. */
	private final TaxIdNode excluded = new TaxIdNode("2608887");
	private final TaxIdNode belowExcluded = new TaxIdNode("1306");

	private Set<TaxIdNode> exclusion() {
		// What ExcludedTaxNodesGoal hands the update: the struck-out tax ids *with* their descendants.
		Set<TaxIdNode> nodes = new HashSet<>();
		nodes.add(excluded);
		nodes.add(belowExcluded);
		return nodes;
	}

	private Set<TaxIdNode> selection() {
		// What TaxNodesGoal hands the update: the requested tax ids *with* their descendants.
		Set<TaxIdNode> taxNodes = new HashSet<>();
		taxNodes.add(requested);
		taxNodes.add(strainBelowRequested);
		return taxNodes;
	}

	/** {@link DBGoal#isRegionInScope} for a region of the release, which is what the scope restricts. */
	private boolean inScope(UpdateScope scope, Set<TaxIdNode> taxNodes, TaxIdNode node) {
		return DBGoal.isRegionInScope(scope, taxNodes, exclusion(), node, FROM_RELEASE);
	}

	@Test
	public void allIncludesEverything() {
		Set<TaxIdNode> taxNodes = selection();
		assertTrue(inScope(UpdateScope.ALL, taxNodes, requested));
		assertTrue(inScope(UpdateScope.ALL, taxNodes, strainBelowRequested));
		assertTrue(inScope(UpdateScope.ALL, taxNodes, other));
		assertTrue(inScope(UpdateScope.ALL, taxNodes, null));
	}

	@Test
	public void ownTaxaOnlyKeepsTheSelectionAndItsDescendants() {
		Set<TaxIdNode> taxNodes = selection();
		assertTrue(inScope(UpdateScope.OWN_TAXA_ONLY, taxNodes, requested));
		assertTrue(inScope(UpdateScope.OWN_TAXA_ONLY, taxNodes, strainBelowRequested));
		assertFalse(inScope(UpdateScope.OWN_TAXA_ONLY, taxNodes, other));
		// A region whose accession is not in the map has no taxon, so it is not one of ours.
		assertFalse(inScope(UpdateScope.OWN_TAXA_ONLY, taxNodes, null));
	}

	@Test
	public void otherTaxaOnlyIsTheComplementOfOwnTaxaOnly() {
		Set<TaxIdNode> taxNodes = selection();
		// This is the point of the scope: the genomes of the requested taxon are already in the
		// database under their own identity, so the release's copy of them must not be merged in.
		assertFalse(inScope(UpdateScope.OTHER_TAXA_ONLY, taxNodes, requested));
		assertFalse(inScope(UpdateScope.OTHER_TAXA_ONLY, taxNodes, strainBelowRequested));
		// And this is what it keeps - without it the update would have no purpose at all.
		assertTrue(inScope(UpdateScope.OTHER_TAXA_ONLY, taxNodes, other));
		assertTrue(inScope(UpdateScope.OTHER_TAXA_ONLY, taxNodes, null));

		for (TaxIdNode node : new TaxIdNode[] { requested, strainBelowRequested, other, null }) {
			assertEquals("complementary for " + node,
					inScope(UpdateScope.OWN_TAXA_ONLY, taxNodes, node),
					!inScope(UpdateScope.OTHER_TAXA_ONLY, taxNodes, node));
		}
	}

	@Test
	public void anEmptySelectionRestrictsNothingInEitherDirection() {
		Set<TaxIdNode> none = Collections.emptySet();
		// An empty selection is "no restriction", and both restricted scopes have to read it that way:
		// ownTaxaOnly by its explicit isEmpty() case, otherTaxaOnly because nothing is then ours.
		assertTrue(inScope(UpdateScope.OWN_TAXA_ONLY, none, requested));
		assertTrue(inScope(UpdateScope.OWN_TAXA_ONLY, none, other));
		assertTrue(inScope(UpdateScope.OTHER_TAXA_ONLY, none, requested));
		assertTrue(inScope(UpdateScope.OTHER_TAXA_ONLY, none, other));
	}

	@Test
	public void projectFastasAreNeverOutOfScope() {
		Set<TaxIdNode> taxNodes = selection();
		// The scope is named refseq.updateScope because it restricts the release and nothing else. A
		// genome the project supplies itself - an additional.txt entry or a Genbank download - is in
		// the database under the very identity the update meets it with, so LCA(n, n) = n leaves it
		// where it is; what the update does do there is settle a k-mer that several of those genomes
		// share on their common ancestor instead of leaving it claimed for whichever was read first.
		// That is the half of the update a per-assembly database lives on, and no scope may drop it.
		for (UpdateScope scope : UpdateScope.values()) {
			for (TaxIdNode node : new TaxIdNode[] { requested, strainBelowRequested, other, null }) {
				assertTrue(scope.getName() + " / " + node,
						DBGoal.isRegionInScope(scope, taxNodes, exclusion(), node, FROM_PROJECT_FASTA));
			}
		}
	}

	@Test
	public void everyScopeIsHandled() {
		Set<TaxIdNode> taxNodes = selection();
		for (UpdateScope scope : UpdateScope.values()) {
			// No scope may fall through to an exception; the default branch has to cover ALL only.
			inScope(scope, taxNodes, requested);
		}
	}

	@Test
	public void theConfigKeyParsesItsNames() {
		ConfigParamInfo<?> info = GSConfigKey.UPDATE_SCOPE.getInfo();
		assertEquals(UpdateScope.ALL, info.defaultValue());
		// The `refseq.' prefix is part of the contract: it says the scope restricts the release only.
		assertEquals("refseq.updateScope", GSConfigKey.UPDATE_SCOPE.getName());
		assertFalse(GSConfigKey.UPDATE_SCOPE.isInternal());

		for (UpdateScope scope : UpdateScope.values()) {
			assertEquals(scope, UpdateScope.byName(scope.getName()));
			assertTrue(scope.getName(), info.isValid(scope.getName()));
			assertTrue(scope.getName(), info.isInRange(scope.getName()));
		}
		// Case is not significant, the constant's own name is not a configuration name, and an
		// unknown value must be rejected rather than silently taken for one of the three.
		assertEquals(UpdateScope.OTHER_TAXA_ONLY, UpdateScope.byName("OTHERTAXAONLY"));
		assertNull(UpdateScope.byName("OTHER_TAXA_ONLY"));
		assertNull(UpdateScope.byName("minUpdate"));
		assertFalse(info.isValid("true"));
	}

	@Test
	public void theDocumentedDefaultIsAValueThatCanBeConfigured() {
		// The generated ConfigParams.md must not name a default that a properties file would reject.
		ConfigParamInfo<?> info = GSConfigKey.UPDATE_SCOPE.getInfo();
		assertEquals("all", info.getMDDefaultValue());
		assertTrue(info.isValid(info.getMDDefaultValue()));
	}

	/**
	 * {@code allButExcluded} is {@code all} with the struck-out branches taken out of it, which is
	 * what {@code all} does not do: an exclusion in {@code taxids.txt} keeps a branch out of the
	 * database but, under {@code all}, its genomes still meet the kept taxa in the update and raise
	 * their k-mers to a common ancestor.
	 */
	@Test
	public void allButExcludedSkipsTheStruckOutBranch() {
		Set<TaxIdNode> taxNodes = selection();
		assertFalse(inScope(UpdateScope.ALL_BUT_EXCLUDED, taxNodes, excluded));
		assertFalse(inScope(UpdateScope.ALL_BUT_EXCLUDED, taxNodes, belowExcluded));
	}

	/** Everything else takes part, which is what separates it from {@code ownTaxaOnly}. */
	@Test
	public void allButExcludedKeepsEverythingElse() {
		Set<TaxIdNode> taxNodes = selection();
		assertTrue(inScope(UpdateScope.ALL_BUT_EXCLUDED, taxNodes, requested));
		assertTrue(inScope(UpdateScope.ALL_BUT_EXCLUDED, taxNodes, strainBelowRequested));
		// The point of the scope: a genome of another taxon still raises k-mers, so a species does
		// not keep a specificity that its relatives disprove.
		assertTrue(inScope(UpdateScope.ALL_BUT_EXCLUDED, taxNodes, other));
	}

	/**
	 * A region whose accession is not in the map has no taxon and cannot be one of the excluded ones,
	 * so it takes part as it does under {@code all}. It cannot raise anything either way.
	 */
	@Test
	public void allButExcludedKeepsARegionOfUnknownTaxon() {
		assertTrue(inScope(UpdateScope.ALL_BUT_EXCLUDED, selection(), null));
	}

	/** With nothing struck out the scope is {@code all}, node for node. */
	@Test
	public void allButExcludedWithoutExclusionsIsAll() {
		Set<TaxIdNode> taxNodes = selection();
		for (TaxIdNode node : new TaxIdNode[] { requested, strainBelowRequested, other, excluded, null }) {
			assertEquals("node " + node,
					DBGoal.isRegionInScope(UpdateScope.ALL, taxNodes, new HashSet<TaxIdNode>(), node, FROM_RELEASE),
					DBGoal.isRegionInScope(UpdateScope.ALL_BUT_EXCLUDED, taxNodes, new HashSet<TaxIdNode>(), node,
							FROM_RELEASE));
		}
	}
}
