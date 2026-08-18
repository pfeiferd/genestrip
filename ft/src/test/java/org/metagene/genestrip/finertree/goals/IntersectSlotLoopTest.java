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

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;

import java.util.Random;

import org.junit.Test;

/**
 * Pins down the loop transformation in {@link KMerIntersectCountGoal#inKMerStoreWork}: counting the
 * pairs among the <em>set</em> slots has to give the same triangular count array as the previous
 * sweep over every pair of slots.
 * <p>
 * The counts this array holds feed the Jaccard indices, hence the dendrogram, hence which k-mer ends
 * up on which node - so "faster" is worth nothing here unless it is also identical. Both loops are
 * written out below and compared over random membership vectors, including the shapes that are easy
 * to get wrong: nothing set, everything set, and only the trailing OTHER slot set.
 * <p>
 * What made the sweep worth replacing is that its cost does not depend on how many slots are set. A
 * node with 3,500 children ran six million iterations per k-mer to count the two or three pairs a
 * typical k-mer contributes.
 */
public class IntersectSlotLoopTest {
    /** The count array as {@code IntersectionsPerNodeImpl} sizes it for a node with that many children. */
    private static long[] newCounts(int children) {
        int c = children + 1;
        return new long[(c * c + c) / 2 + 2];
    }

    /** The sweep over every pair of slots, with the index normalisation the old path applied. */
    private static long[] byFullSweep(boolean[] bits, int children, int spread) {
        long[] counts = newCounts(children);
        counts[counts.length - 2] += spread;
        counts[counts.length - 1]++;
        for (int i = 0; i <= children; i++) {
            for (int j = i; j <= children; j++) {
                if (bits[i] && bits[j]) {
                    int lo = Math.min(i, j);
                    int hi = Math.max(i, j);
                    counts[(hi * hi + hi) / 2 + lo]++;
                }
            }
        }
        return counts;
    }

    /** The loop over the set slots alone, as the goal runs it now. */
    private static long[] bySetSlots(int[] setSlots, int children, int spread) {
        long[] counts = newCounts(children);
        counts[counts.length - 2] += spread;
        counts[counts.length - 1]++;
        for (int a = 0; a < spread; a++) {
            int i = setSlots[a];
            for (int b = a; b < spread; b++) {
                int j = setSlots[b];
                counts[(j * j + j) / 2 + i]++;
            }
        }
        return counts;
    }

    /** Fills the slot list the way the visitor's membership sweep does, ascending. */
    private static int collect(boolean[] bits, int children, int[] setSlots) {
        int spread = 0;
        for (int i = 0; i <= children; i++) {
            if (bits[i]) {
                setSlots[spread++] = i;
            }
        }
        return spread;
    }

    private void assertSame(boolean[] bits, int children) {
        int[] setSlots = new int[children + 1];
        int spread = collect(bits, children, setSlots);
        assertArrayEquals("children=" + children + " spread=" + spread,
                byFullSweep(bits, children, spread), bySetSlots(setSlots, children, spread));
    }

    @Test
    public void testTheTwoLoopsAgreeOnRandomMembership() {
        // A fixed seed: a test that fails only sometimes says less than one that fails always.
        Random random = new Random(20260814L);
        for (int children = 1; children <= 40; children++) {
            for (int round = 0; round < 50; round++) {
                boolean[] bits = new boolean[children + 1];
                // Sparse on purpose - that is the case the change exists for, a handful of the
                // children carrying the k-mer.
                for (int i = 0; i <= children; i++) {
                    bits[i] = random.nextInt(10) == 0;
                }
                assertSame(bits, children);
            }
        }
    }

    @Test
    public void testNothingSetCountsNoPair() {
        boolean[] bits = new boolean[9];
        assertSame(bits, 8);
        int[] setSlots = new int[9];
        int spread = collect(bits, 8, setSlots);
        assertEquals(0, spread);
        long[] counts = bySetSlots(setSlots, 8, spread);
        for (int i = 0; i < counts.length - 2; i++) {
            assertEquals("no pair may be counted", 0, counts[i]);
        }
        assertEquals("but the k-mer is still one k-mer", 1, counts[counts.length - 1]);
    }

    @Test
    public void testEverythingSetIsTheWorstCaseAndStillAgrees() {
        boolean[] bits = new boolean[13];
        java.util.Arrays.fill(bits, true);
        assertSame(bits, 12);
    }

    /**
     * Only the trailing OTHER slot set. It is the highest index and the one a stray {@code <=} would
     * drop, and it means the k-mer is carried by something that is none of the children.
     */
    @Test
    public void testOnlyTheOtherSlotSet() {
        boolean[] bits = new boolean[11];
        bits[10] = true;
        assertSame(bits, 10);
        int[] setSlots = new int[11];
        int spread = collect(bits, 10, setSlots);
        assertEquals(1, spread);
        assertEquals(10, setSlots[0]);
        long[] counts = bySetSlots(setSlots, 10, spread);
        assertEquals("the OTHER slot paired with itself", 1, counts[(10 * 10 + 10) / 2 + 10]);
    }

    /** A single child and the OTHER slot - the smallest node the refinement ever considers. */
    @Test
    public void testTheSmallestNode() {
        boolean[] bits = { true, true };
        assertSame(bits, 1);
        boolean[] onlyChild = { true, false };
        assertSame(onlyChild, 1);
    }
}
