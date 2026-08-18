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
 * The first is the OTHER slot. It sits one past the children, says the k-mer was read from a genome
 * below none of them, and no refined node covers it - so such a k-mer has to stay where it is. The
 * bit sets are only as long as the child list, so a comparison bounded by the <em>container</em>
 * never looks at that slot and pushes the k-mer down anyway. That was the behaviour until
 * 2026-08-14, and it made the refinement claim a specificity the genome in the OTHER slot
 * contradicts. {@link #testTheOtherSlotBlocksEveryRefinedNode()} is the regression.
 * <p>
 * The second is that {@code bits} is the visitor's buffer, sized to a power of two and reused for
 * every parent, so its tail holds whatever a previous, possibly larger, parent left there. The
 * length has to come from the parent being refined and not from the array -
 * {@link #testStaleSlotsBeyondThisParentAreNotRead()}.
 */
public class UpdateStoreCoversAllTest {
    /** Slots of a parent with four children: four children plus OTHER. */
    private static final int LEN = 5;

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
        assertTrue(UpdateStoreGoal.coversAll(bits("1100"), bits("10000"), LEN));
        assertTrue(UpdateStoreGoal.coversAll(bits("1100"), bits("11000"), LEN));
        assertTrue("covering more than needed is still covering",
                UpdateStoreGoal.coversAll(bits("1111"), bits("01000"), LEN));
    }

    /** A node missing one of them does not. */
    @Test
    public void testANodeMissingAChildDoesNot() {
        assertFalse(UpdateStoreGoal.coversAll(bits("1100"), bits("10100"), LEN));
        assertFalse(UpdateStoreGoal.coversAll(bits("0011"), bits("11000"), LEN));
    }

    /**
     * The regression. With the OTHER slot set, no refined node may take the k-mer - not even one that
     * covers every child, and not even the one covering all of them.
     */
    @Test
    public void testTheOtherSlotBlocksEveryRefinedNode() {
        assertFalse("a node covering all four children still must not take it",
                UpdateStoreGoal.coversAll(bits("1111"), bits("11111"), LEN));
        assertFalse(UpdateStoreGoal.coversAll(bits("1100"), bits("10001"), LEN));
        assertFalse("OTHER alone is enough to block it",
                UpdateStoreGoal.coversAll(bits("1111"), bits("00001"), LEN));
    }

    /**
     * The buffer is reused across parents, so slots past this parent's own must not be read whatever
     * they happen to hold.
     */
    @Test
    public void testStaleSlotsBeyondThisParentAreNotRead() {
        boolean[] reused = bits("10000" + "1111111");   // 5 live slots, 7 slots of another parent's
        assertTrue(UpdateStoreGoal.coversAll(bits("1100"), reused, LEN));
        boolean[] staleOther = bits("00000" + "1");
        assertTrue("a stale bit must not be mistaken for OTHER",
                UpdateStoreGoal.coversAll(bits("1000"), staleOther, LEN));
    }

    /** A k-mer found in nothing at all is covered by anything; the caller never asks, but it is total. */
    @Test
    public void testNothingSetIsCoveredByAnything() {
        assertTrue(UpdateStoreGoal.coversAll(bits("0000"), bits("00000"), LEN));
    }

    /**
     * And the selection itself: the bit sets are kept sorted by ascending cardinality, so the first
     * covering one is the most specific. This walks them the way {@code getBestMatchingNode} does and
     * checks that the answer really is the smallest node that covers the k-mer.
     */
    @Test
    public void testTheFirstCoveringSetIsTheMostSpecificOne() {
        // Sorted by cardinality, as BitSetsForNodes.sort() leaves them.
        boolean[][] sets = { bits("1100"), bits("0011"), bits("1110"), bits("1111") };
        assertEquals("the pair, not the triple or the root", 0, firstCovering(sets, bits("10000")));
        assertEquals(0, firstCovering(sets, bits("11000")));
        assertEquals("needs a child of each pair, so the triple", 2, firstCovering(sets, bits("10100")));
        assertEquals("spans both pairs, so the root", 3, firstCovering(sets, bits("10010")));
        assertEquals("OTHER set: nothing takes it", -1, firstCovering(sets, bits("10001")));
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
