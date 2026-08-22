package org.metagene.genestrip.finertree.goals;

import org.junit.Test;
import org.metagene.genestrip.finertree.probfilter.KMerIndexFilterHelper;
import org.metagene.genestrip.tax.Rank;
import org.metagene.genestrip.tax.SmallTaxTree.SmallTaxIdNode;

import java.util.HashSet;
import java.util.Set;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotEquals;
import static org.junit.Assert.assertTrue;

/**
 * Pins the seam the two goal families share and disagree on: which key a (<em>k</em>-mer, node) pair
 * folds into.
 * <p>
 * Both read the same sequences, resolve the same leaf and combine it with the same {@code k}-mer
 * through {@link KMerIndexFilterHelper#combine}. What they put in as the second component differs,
 * and that difference is the whole of what separates them:
 * <ul>
 * <li>the index family ({@code kmerindexsize}, {@code kmerindexbloom}) uses the <em>index of the
 *     direct child of the refined node</em> the leaf lies under, or {@code OTHER_VALUE};</li>
 * <li>the quality family ({@code dbqualcounts} and its sizing pass) uses the <em>leaf's own
 *     position</em>.</li>
 * </ul>
 * Each choice is right for its purpose and wrong for the other's. The index asks which child subtree
 * carries a {@code k}-mer, so two leaves under one child must collapse to one key; the quality
 * measure asks how many data taxa carry it, so those same two leaves must stay apart. Swapping them
 * would not fail loudly: it would report a database as better or worse placed than it is, which is
 * exactly the kind of defect these measurements exist to detect elsewhere.
 * <p>
 * {@link KMerIndexOtherSlotTest} pins the child index and {@link
 * org.metagene.genestrip.finertree.probfilter.KMerIndexFilterHelperTest} the folding. This pins that
 * the two are composed the way each family needs, which no other test does.
 */
public class PairKeyCompositionTest {
    private final SmallTaxIdNode file1 = leaf("file1", 11);
    private final SmallTaxIdNode file2 = leaf("file2", 12);
    private final SmallTaxIdNode file3 = leaf("file3", 13);
    private final SmallTaxIdNode data = inner("data", 2, Rank.DATA, file1, file2);
    private final SmallTaxIdNode data2 = inner("data2", 4, Rank.DATA, file3);
    private final SmallTaxIdNode strain = inner("strain", 3, Rank.STRAIN, data2);
    private final SmallTaxIdNode species = inner("species", 1, Rank.SPECIES, data, strain);
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

    /** What the index family puts into the filter for one pair. */
    private static long indexKey(long kmer, SmallTaxIdNode leaf, SmallTaxIdNode storedNode) {
        return KMerIndexFilterHelper.combine(kmer, AbstractKMerIndexGoal.childIndexUnder(leaf, storedNode));
    }

    /**
     * What the quality family puts into the filter for one pair. The position is passed in rather
     * than read off the node: a {@link SmallTaxIdNode} gets its dense position when a whole
     * {@link org.metagene.genestrip.tax.SmallTaxTree} is built, and the nodes here are hand-made.
     * What is pinned is that the quality family keys by the leaf and the index family by the child.
     */
    private static long qualityKey(long kmer, int leafPos) {
        return KMerIndexFilterHelper.combine(kmer, leafPos);
    }

    @Test
    public void testTwoLeavesUnderOneChildCollapseForTheIndex() {
        // file1 and file2 both hang under `data', which is the direct child of `species'. The index
        // asks which child subtree carries the k-mer, so one k-mer in both is one entry.
        long kmer = 0x0123456789ABCDEFL;
        assertEquals("the index counts child subtrees, not leaves",
                indexKey(kmer, file1, species), indexKey(kmer, file2, species));
    }

    @Test
    public void testTheSameTwoLeavesStayApartForTheQualityMeasure() {
        // The quality measure asks how many data taxa carry the k-mer, so the two must not collapse,
        // or a k-mer in two genomes would count as one and every precision derived from it would be
        // wrong in the direction of looking better.
        long kmer = 0x0123456789ABCDEFL;
        int file1Pos = 7;
        int file2Pos = 8;
        assertNotEquals("the quality measure counts leaves, not child subtrees",
                qualityKey(kmer, file1Pos), qualityKey(kmer, file2Pos));
    }

    @Test
    public void testADeeperLeafFoldsOntoItsChildForTheIndexButKeepsItselfForQuality() {
        long kmer = 42L;
        int file3Pos = 21;
        assertEquals("file3 lies under `strain', the child of `species'",
                KMerIndexFilterHelper.combine(kmer, strain.storeIndex), indexKey(kmer, file3, species));
        assertEquals("the quality measure takes the leaf itself",
                KMerIndexFilterHelper.combine(kmer, file3Pos), qualityKey(kmer, file3Pos));
    }

    @Test
    public void testALeafOutsideTheNodeGoesToTheOtherSlotForTheIndex() {
        // A k-mer the LCA update pushed up because of a taxon the database does not hold: the index
        // needs the OTHER slot for it, or the refinement would push the k-mer down to the one species
        // it can see. The quality measure has no such case - it only ever sees leaves of its own tree.
        long kmer = 99L;
        assertEquals(KMerIndexFilterHelper.combine(kmer, AbstractKMerIndexGoal.OTHER_VALUE),
                indexKey(kmer, elsewhere, species));
    }

    @Test
    public void testKeysSeparateBothComponents() {
        // The folding must keep two k-mers apart under one node and one k-mer apart under two nodes.
        // combine() is an xor with a multiple of the index, so neither is self-evident: an index
        // whose product collided with a k-mer difference would silently merge two pairs into one and
        // undersize every filter built from the count.
        Set<Long> keys = new HashSet<>();
        int collisions = 0;
        for (long kmer = 0; kmer < 500; kmer++) {
            for (int pos = 0; pos < 500; pos++) {
                if (!keys.add(KMerIndexFilterHelper.combine(kmer, pos))) {
                    collisions++;
                }
            }
        }
        assertEquals("250,000 distinct pairs must fold to 250,000 distinct keys", 0, collisions);
        assertEquals(250_000, keys.size());
    }

    @Test
    public void testTheOtherSlotDoesNotCollideWithARealIndex() {
        // OTHER_VALUE is Integer.MAX_VALUE, i.e. a store index no node will ever have - but the keys
        // are what matter, and they are products. A collision here would attribute k-mers of an
        // organism outside the database to a species inside it.
        long kmer = 0x5555555555555555L;
        long other = KMerIndexFilterHelper.combine(kmer, AbstractKMerIndexGoal.OTHER_VALUE);
        for (int pos = 0; pos < 10_000; pos++) {
            assertNotEquals("store index " + pos + " collides with the OTHER slot",
                    other, KMerIndexFilterHelper.combine(kmer, pos));
        }
    }

    @Test
    public void testTheTwoFamiliesDisagreeWhereTheyMust() {
        // The point of the whole test: for a leaf that is not itself the direct child, the two keys
        // differ. A refactor that gave one family the other's second component would pass every
        // existing test and change every number.
        long kmer = 7L;
        int file3Pos = 21;
        assertNotEquals("child index and leaf position must not be interchangeable here",
                indexKey(kmer, file3, species), qualityKey(kmer, file3Pos));
        assertTrue("file3 sits below the child, not at it", strain.storeIndex != file3Pos);
    }
}
