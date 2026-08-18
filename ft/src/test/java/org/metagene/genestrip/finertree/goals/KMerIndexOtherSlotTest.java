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
package org.metagene.genestrip.finertree.goals;

import static org.junit.Assert.assertEquals;

import org.junit.Test;
import org.metagene.genestrip.tax.Rank;
import org.metagene.genestrip.tax.SmallTaxTree.SmallTaxIdNode;

/**
 * Tests which slot of the k-mer index a (k-mer, leaf) pair is recorded under.
 * <p>
 * The reassignment asks the filter once per <em>direct</em> child of the node being refined and once
 * for the trailing OTHER slot, and for nothing else. So that is what the index has to hold: the
 * genome a k-mer was read from is mapped onto the child it lies under, however deep below that child
 * it sits. Recording the genome itself would fill the index with entries nothing ever asks for - one
 * per genome instead of one per child - and would take the bound {@code kmerindexbloom} sizes its
 * filter by with it, since that bound counts children.
 * <p>
 * The OTHER slot catches everything that is not one of the children: a genome that resolved to the
 * refined node itself, and one that is not below it at all. A pair recorded under the node's own
 * store index instead would be queried by nothing, and the k-mer would look as though that genome
 * did not carry it - the very thing OTHER exists to prevent.
 */
public class KMerIndexOtherSlotTest {
	/**
	 * The shape a per-assembly database has below the node it refines:
	 *
	 * <pre>
	 * species -&gt; data -&gt; file1, file2
	 *         -&gt; strain -&gt; data2 -&gt; file3
	 * </pre>
	 */
	private final SmallTaxIdNode file1 = leaf("file1", 11);
	private final SmallTaxIdNode file2 = leaf("file2", 12);
	private final SmallTaxIdNode file3 = leaf("file3", 13);
	private final SmallTaxIdNode data = inner("data", 2, Rank.DATA, file1, file2);
	private final SmallTaxIdNode data2 = inner("data2", 4, Rank.DATA, file3);
	private final SmallTaxIdNode strain = inner("strain", 3, Rank.STRAIN, data2);
	private final SmallTaxIdNode species = inner("species", 1, Rank.SPECIES, data, strain);
	/** A node of another lineage, so not below {@link #species} at all. */
	private final SmallTaxIdNode elsewhere = leaf("elsewhere", 99);

	private static SmallTaxIdNode leaf(String taxId, int storeIndex) {
		SmallTaxIdNode node = new SmallTaxIdNode(taxId, taxId, Rank.FILE);
		node.setStoreIndex(storeIndex);
		return node;
	}

	private static SmallTaxIdNode inner(String taxId, int storeIndex, Rank rank, SmallTaxIdNode... subNodes) {
		SmallTaxIdNode node = new SmallTaxIdNode(taxId, taxId, rank, subNodes);
		node.setStoreIndex(storeIndex);
		return node;
	}

	/** A direct child is recorded as itself; that is the ordinary case. */
	@Test
	public void testADirectChildKeepsItsIndex() {
		assertEquals(data.storeIndex, AbstractKMerIndexGoal.childIndexUnder(data, species));
		assertEquals(strain.storeIndex, AbstractKMerIndexGoal.childIndexUnder(strain, species));
		assertEquals(file1.storeIndex, AbstractKMerIndexGoal.childIndexUnder(file1, data));
	}

	/** A leaf further down is recorded under the child it lies below, not under itself. */
	@Test
	public void testADeeperLeafIsMappedOntoTheChildItLiesUnder() {
		// This is the whole point: three assemblies, one entry each under `data' or `strain', and
		// never one under file1, file2 or file3 - which nothing would ever query.
		assertEquals(data.storeIndex, AbstractKMerIndexGoal.childIndexUnder(file1, species));
		assertEquals(data.storeIndex, AbstractKMerIndexGoal.childIndexUnder(file2, species));
		assertEquals(strain.storeIndex, AbstractKMerIndexGoal.childIndexUnder(file3, species));
		// Two levels further down still resolves to the same direct child.
		assertEquals(strain.storeIndex, AbstractKMerIndexGoal.childIndexUnder(data2, species));
	}

	/** A leaf that is the refined node itself belongs in the OTHER slot. */
	@Test
	public void testTheNodeItselfGoesToTheOtherSlot() {
		assertEquals(AbstractKMerIndexGoal.OTHER_VALUE, AbstractKMerIndexGoal.childIndexUnder(species, species));
		assertEquals(AbstractKMerIndexGoal.OTHER_VALUE, AbstractKMerIndexGoal.childIndexUnder(file1, file1));
	}

	/**
	 * A leaf outside the refined node's subtree goes to OTHER as well. The lowest common ancestor
	 * rules that case out, but a scope that withholds genomes from the update does not - and then
	 * OTHER is both the true answer and what keeps the index bounded by the children.
	 */
	@Test
	public void testALeafOutsideTheSubtreeGoesToTheOtherSlot() {
		assertEquals(AbstractKMerIndexGoal.OTHER_VALUE, AbstractKMerIndexGoal.childIndexUnder(elsewhere, species));
		// Upwards is outside too: the species is not below its own data node.
		assertEquals(AbstractKMerIndexGoal.OTHER_VALUE, AbstractKMerIndexGoal.childIndexUnder(species, data));
	}

	/** A region that resolved to no leaf at all is recorded under OTHER. */
	@Test
	public void testNoLeafGoesToTheOtherSlot() {
		assertEquals(AbstractKMerIndexGoal.OTHER_VALUE, AbstractKMerIndexGoal.childIndexUnder(null, species));
	}
}
