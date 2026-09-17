/*
 * Genestrip
 */
package org.metagene.genestrip.finertree.goals;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.junit.Test;
import org.metagene.genestrip.finertree.cluster.DendrogramNode;

/**
 * Tests the bit sets the refinement builds from a dendrogram, against dendrograms rather than against
 * hand-written vectors.
 * <p>
 * This exists because {@link UpdateStoreCoversAllTest} could not catch what went wrong here. That
 * test exercises the comparison on vectors it writes itself, and the comparison was right the whole
 * time; what was wrong was the shape of the vectors the refinement built - one slot per child and
 * none for OTHER, so that a refined node whose cluster contained OTHER could never be seen to cover
 * it. Every k-mer an unnamed genome touches was stranded at the parent, and the unit test went on
 * passing. A test of a comparison cannot find a defect in its inputs, so this one starts from the
 * dendrogram.
 * <p>
 * The measured case: on the Orthopoxvirus database of the companion paper the genus holds 239,260
 * k-mers, of which 97,474 are touched by a genome the database does not name. Those stayed at the
 * genus, which came out at 233,212 where the published figure is 134,449.
 */
public class UpdateStoreBitSetsTest {

    /** A leaf standing for the child at that index; the OTHER leaf is the one at {@code childCount}. */
    private static DendrogramNode leaf(int valueIndex) {
        return new DendrogramNode(valueIndex, 1.0);
    }

    private static DendrogramNode merge(DendrogramNode a, DendrogramNode b) {
        return new DendrogramNode(a, b, 0.5);
    }

    private static boolean[] bits(String pattern) {
        boolean[] b = new boolean[pattern.length()];
        for (int i = 0; i < pattern.length(); i++) {
            b[i] = pattern.charAt(i) == '1';
        }
        return b;
    }

    /**
     * Every bit set has a slot per child and one for OTHER. This is the defect itself: the sets used
     * to be one slot short, and everything below follows from that.
     */
    @Test
    public void testABitSetHasASlotForOther() {
        // Three children V, C, M and the OTHER leaf at index 3. Clustered ((V, OTHER), C) beside M,
        // which is what the Orthopoxvirus genus actually produced.
        DendrogramNode vOther = merge(leaf(0), leaf(3));
        DendrogramNode root = merge(merge(vOther, leaf(1)), leaf(2));
        for (boolean[] set : UpdateStoreGoal.buildBitSets(3, root)) {
            assertEquals("a bit set must have one slot per child plus OTHER", 4, set.length);
        }
    }

    /**
     * A refined node whose cluster contains OTHER records it, and one whose cluster does not leaves it
     * clear. Without this the two cases are indistinguishable and the refinement can only refuse both.
     */
    @Test
    public void testTheOtherLeafSetsItsSlot() {
        DendrogramNode vOther = merge(leaf(0), leaf(3));
        DendrogramNode root = merge(merge(vOther, leaf(1)), leaf(2));
        boolean[][] sets = UpdateStoreGoal.buildBitSets(3, root);
        assertTrue("{V, OTHER} must be among the sets built", contains(sets, bits("1001")));
        assertTrue("{V, C, OTHER} must be too", contains(sets, bits("1101")));

        // The same shape with two real children in place of V and OTHER: nothing sets the last slot.
        DendrogramNode noOther = merge(merge(merge(leaf(0), leaf(1)), leaf(2)), leaf(3));
        for (boolean[] set : UpdateStoreGoal.buildBitSets(4, noOther)) {
            assertFalse("no cluster contains OTHER, so no set may claim it", set[4]);
        }
    }

    /**
     * And the consequence, which is what the database gets wrong or right: a k-mer carried by a child
     * and by an unnamed genome is taken by the node clustering that child with OTHER, and a k-mer
     * carried by a child the node does not cover is not.
     */
    @Test
    public void testAnOtherKMerReachesTheNodeThatCoversIt() {
        DendrogramNode vOther = merge(leaf(0), leaf(3));
        DendrogramNode root = merge(merge(vOther, leaf(1)), leaf(2));
        boolean[][] sets = UpdateStoreGoal.buildBitSets(3, root);

        assertTrue("carried by V and an unnamed genome: {V, OTHER} takes it",
                covers(sets, bits("1001"), bits("1001")));
        assertTrue("carried by V and C and an unnamed genome: {V, C, OTHER} takes it",
                covers(sets, bits("1101"), bits("1101")));
        assertTrue("carried by V and C alone: {V, C, OTHER} covers it as well",
                covers(sets, bits("1101"), bits("1100")));
        assertFalse("carried by M as well: no refined node covers M, so it stays at the genus",
                anyCovers(sets, bits("1011")));
    }

    /**
     * Where no cluster contains OTHER, a k-mer an unnamed genome touches has nowhere to go and must
     * stay at the parent. That was the original defect, fixed on 2026-08-14, and it must not come
     * back while the opposite one is being fixed.
     */
    @Test
    public void testWithoutSuchANodeAnOtherKMerStaysAtTheParent() {
        DendrogramNode noOther = merge(merge(merge(leaf(0), leaf(1)), leaf(2)), leaf(3));
        boolean[][] sets = UpdateStoreGoal.buildBitSets(4, noOther);
        assertTrue("without the OTHER bit it is taken", anyCovers(sets, bits("11000")));
        assertFalse("with it, nothing takes it", anyCovers(sets, bits("11001")));
        assertFalse("not even a k-mer only an unnamed genome carries", anyCovers(sets, bits("00001")));
    }

    /** The union is what a set records: a parent covers exactly what its two children cover. */
    @Test
    public void testASetIsTheUnionOfItsChildren() {
        DendrogramNode inner = merge(leaf(0), leaf(3));
        DendrogramNode root = merge(merge(inner, leaf(1)), leaf(2));
        boolean[][] sets = UpdateStoreGoal.buildBitSets(3, root);
        boolean[] vOther = find(sets, 2);
        assertArrayEquals("V and OTHER and nothing else", bits("1001"), vOther);
    }

    private static boolean contains(boolean[][] sets, boolean[] wanted) {
        for (boolean[] set : sets) {
            if (java.util.Arrays.equals(set, wanted)) {
                return true;
            }
        }
        return false;
    }

    /** The set equal to {@code wanted} must cover {@code kmer}. */
    private static boolean covers(boolean[][] sets, boolean[] wanted, boolean[] kmer) {
        for (boolean[] set : sets) {
            if (java.util.Arrays.equals(set, wanted)) {
                return UpdateStoreGoal.coversAll(set, kmer, kmer.length);
            }
        }
        return false;
    }

    private static boolean anyCovers(boolean[][] sets, boolean[] kmer) {
        for (boolean[] set : sets) {
            if (UpdateStoreGoal.coversAll(set, kmer, kmer.length)) {
                return true;
            }
        }
        return false;
    }

    /** The first set of the given cardinality. */
    private static boolean[] find(boolean[][] sets, int cardinality) {
        for (boolean[] set : sets) {
            int c = 0;
            for (boolean b : set) {
                if (b) {
                    c++;
                }
            }
            if (c == cardinality) {
                return set;
            }
        }
        return null;
    }
}
