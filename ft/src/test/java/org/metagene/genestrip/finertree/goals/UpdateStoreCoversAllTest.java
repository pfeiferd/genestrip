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
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

/**
 * Tests the predicate that decides where a k-mer is pushed to,
 * {@link UpdateStoreGoal#coversAll(boolean[], boolean[], int)}, and the selection built on it.
 * <p>
 * A refined node may take a k-mer only if it covers <em>every</em> slot the k-mer was found in. Two
 * things make that easy to get subtly wrong, and both are pinned down here.
 * <p>
 * The first is the OTHER slot. It sits one past the children and says the k-mer was read from a
 * genome below none of them. It is a slot like any other and both its cases matter. A refined node
 * whose cluster leaves OTHER out does not cover it, and such a k-mer must stay where it is rather
 * than claim a specificity the unnamed genome contradicts -
 * {@link #testANodeWithoutOtherDoesNotTakeAnOtherKMer()}. A refined node whose cluster contains
 * OTHER does cover it, and the k-mer belongs there -
 * {@link #testANodeCoveringOtherTakesAnOtherKMer()}.
 * <p>
 * Both readings have been wrong here before. Until 2026-08-14 the comparison was bounded by the
 * <em>container</em>, which never looked at the slot and pushed every such k-mer down. The fix made
 * OTHER block every refined node instead, which stranded at the parent every k-mer an unnamed genome
 * touches - on the Orthopoxvirus database of the companion paper, ninety-seven thousand of the
 * genus's two hundred and thirty-nine thousand, so that the genus fell to $233{,}212$ where the
 * published figure is $134{,}449$. The bit sets now carry the slot and it is compared like the
 * rest.
 * <p>
 * The second is that {@code bits} is the visitor's buffer, sized to a power of two and reused for
 * every parent, so its tail holds whatever a previous, possibly larger, parent left there. The
 * length has to come from the parent being refined and not from the array -
 * {@link #testStaleSlotsBeyondThisParentAreNotRead()}.
 */
public class UpdateStoreCoversAllTest {
    /** Slots of a parent with four children: four children plus OTHER. */
    private static final int LEN = 5;

    /**
     * Bit sets carry the OTHER slot too, so a container is as long as the vector it is compared
     * against. {@code BitSetsForNodes} sizes them {@code children + 1} for exactly this reason.
     */

    private static boolean[] bits(String pattern) {
        boolean[] b = new boolean[pattern.length()];
        for (int i = 0; i < pattern.length(); i++) {
            b[i] = pattern.charAt(i) == '1';
        }
        return b;
    }

    /** The ordinary case: a node covering the children the k-mer was found in takes it. */
    @Test
    public void testANodeCoveringTheChildrenTakesTheKMer() {
        assertTrue(UpdateStoreGoal.coversAll(bits("11000"), bits("10000"), LEN));
        assertTrue(UpdateStoreGoal.coversAll(bits("11000"), bits("11000"), LEN));
        assertTrue("covering more than needed is still covering",
                UpdateStoreGoal.coversAll(bits("11110"), bits("01000"), LEN));
    }

    /** A node missing one of them does not. */
    @Test
    public void testANodeMissingAChildDoesNot() {
        assertFalse(UpdateStoreGoal.coversAll(bits("11000"), bits("10100"), LEN));
        assertFalse(UpdateStoreGoal.coversAll(bits("00110"), bits("11000"), LEN));
    }

    /**
     * A refined node whose cluster leaves OTHER out must not take a k-mer an unnamed genome touches,
     * however many children it covers. This is the case the parent keeps.
     */
    @Test
    public void testANodeWithoutOtherDoesNotTakeAnOtherKMer() {
        assertFalse("covering all four children is not enough when OTHER is set",
                UpdateStoreGoal.coversAll(bits("11110"), bits("11111"), LEN));
        assertFalse(UpdateStoreGoal.coversAll(bits("11000"), bits("10001"), LEN));
        assertFalse("OTHER alone is enough to block a node that does not cover it",
                UpdateStoreGoal.coversAll(bits("11110"), bits("00001"), LEN));
    }

    /**
     * And the case the blanket refusal destroyed: where the clustering put OTHER inside a refined
     * node, that node covers it and the k-mer belongs there. The node means \"this child or an
     * unnamed relative of it\", which is more specific than the parent and is the whole purpose of
     * computing the slot.
     */
    @Test
    public void testANodeCoveringOtherTakesAnOtherKMer() {
        assertTrue("a child clustered with OTHER takes a k-mer of both",
                UpdateStoreGoal.coversAll(bits("10001"), bits("10001"), LEN));
        assertTrue("and one carried by OTHER alone",
                UpdateStoreGoal.coversAll(bits("10001"), bits("00001"), LEN));
        assertTrue("a larger cluster containing OTHER covers it too",
                UpdateStoreGoal.coversAll(bits("11001"), bits("10001"), LEN));
        assertFalse("but it still needs every child the k-mer was found in",
                UpdateStoreGoal.coversAll(bits("10001"), bits("11001"), LEN));
    }

    /**
     * The buffer is reused across parents, so slots past this parent's own must not be read whatever
     * they happen to hold.
     */
    @Test
    public void testStaleSlotsBeyondThisParentAreNotRead() {
        boolean[] reused = bits("10000" + "1111111");   // 5 live slots, 7 slots of another parent's
        assertTrue(UpdateStoreGoal.coversAll(bits("11000"), reused, LEN));
        boolean[] staleOther = bits("00000" + "1");
        assertTrue("a stale bit must not be mistaken for OTHER",
                UpdateStoreGoal.coversAll(bits("10000"), staleOther, LEN));
    }

    /** A k-mer found in nothing at all is covered by anything; the caller never asks, but it is total. */
    @Test
    public void testNothingSetIsCoveredByAnything() {
        assertTrue(UpdateStoreGoal.coversAll(bits("00000"), bits("00000"), LEN));
    }

    /**
     * And the selection itself: the bit sets are kept sorted by ascending cardinality, so the first
     * covering one is the most specific. This walks them the way {@code getBestMatchingNode} does and
     * checks that the answer really is the smallest node that covers the k-mer.
     */
    @Test
    public void testTheFirstCoveringSetIsTheMostSpecificOne() {
        // Sorted by cardinality, as BitSetsForNodes.sort() leaves them.
        boolean[][] sets = { bits("11000"), bits("00110"), bits("11100"), bits("11110") };
        assertEquals("the pair, not the triple or the root", 0, firstCovering(sets, bits("10000")));
        assertEquals(0, firstCovering(sets, bits("11000")));
        assertEquals("needs a child of each pair, so the triple", 2, firstCovering(sets, bits("10100")));
        assertEquals("spans both pairs, so the root", 3, firstCovering(sets, bits("10010")));
        assertEquals("OTHER set and no set covers it: nothing takes it", -1,
                firstCovering(sets, bits("10001")));
    }

    private static int firstCovering(boolean[][] sets, boolean[] kmerBits) {
        for (int i = 0; i < sets.length; i++) {
            if (UpdateStoreGoal.coversAll(sets[i], kmerBits, LEN)) {
                return i;
            }
        }
        return -1;
    }
}
