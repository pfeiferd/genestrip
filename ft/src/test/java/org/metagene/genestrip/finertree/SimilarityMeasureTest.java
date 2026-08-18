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
import static org.junit.Assert.assertTrue;

import org.junit.Test;
import org.metagene.genestrip.make.ConfigParamInfo;

/**
 * Pins down what the two similarity measures of {@link FTConfigKey#JACCARD_SIM} do to a pair
 * of children, and therefore why one of them is the default and the other a per-database choice.
 * <p>
 * Both are computed from the same three numbers the intersection counts provide - the shared k-mers
 * and each child's own - with Laplace smoothing applied as {@code getJaccardIndex} applies it. The
 * formulas are repeated here rather than driven through the goal, which needs a database.
 */
public class SimilarityMeasureTest {
    private static final double EPS = 1e-9;

    private static final boolean JACCARD = true;
    private static final boolean CONTAINMENT = false;

    private static double similarity(boolean jaccard, long shared, long ownI, long ownJ) {
        long intersect = shared + 1;
        long sizeI = ownI + 1;
        long sizeJ = ownJ + 1;
        return jaccard ? ((double) intersect) / (sizeI + sizeJ - intersect)
                : ((double) intersect) / Math.min(sizeI, sizeJ);
    }

    private static String name(boolean jaccard) {
        return jaccard ? "jaccard" : "containment";
    }

    /** With equally large sets the two agree on identity and differ only in how they fall off. */
    @Test
    public void testTheMeasuresAgreeWhenTheSetsAreTheSameSize() {
        assertEquals(1.0, similarity(JACCARD, 999, 999, 999), EPS);
        assertEquals(1.0, similarity(CONTAINMENT, 999, 999, 999), EPS);
    }

    /**
     * The case containment exists for: two genomes of one lineage, one of them a draft whose k-mer
     * set is a quarter of the other's because the rest was not assembled. Jaccard calls them a
     * quarter similar; containment sees that the draft is contained in the finished genome.
     */
    @Test
    public void testAFragmentedGenomeIsNotPushedAwayByContainment() {
        long finished = 4_000_000;
        long draft = 1_000_000;
        long shared = draft;                       // the draft's k-mers are a subset
        assertEquals(0.25, similarity(JACCARD, shared, finished, draft), 1e-6);
        assertEquals(1.00, similarity(CONTAINMENT, shared, finished, draft), 1e-6);
    }

    /**
     * And the case that keeps Jaccard the default. Where the children are taxa, a small k-mer set is
     * a fact about the taxon and not about how it was sequenced - and containment then scores a taxon
     * with almost nothing stored as maximally similar to whichever one contains it, which would
     * collapse every sparse child onto the densest one.
     */
    @Test
    public void testContainmentCallsATinyChildIdenticalToAHugeOne() {
        long huge = 5_000_000;
        long tiny = 3;
        assertEquals(1.0, similarity(CONTAINMENT, tiny, huge, tiny), EPS);
        assertTrue("Jaccard keeps them far apart",
                similarity(JACCARD, tiny, huge, tiny) < 1e-5);
    }

    /** Sharing nothing is the least similar either measure can report, and neither divides by zero. */
    @Test
    public void testDisjointChildrenAndEmptyChildrenAreSafe() {
        for (boolean jaccard : new boolean[] { JACCARD, CONTAINMENT }) {
            double disjoint = similarity(jaccard, 0, 1000, 1000);
            assertTrue(name(jaccard) + " disjoint", disjoint > 0 && disjoint < 0.01);
            double empty = similarity(jaccard, 0, 0, 0);
            assertEquals(name(jaccard) + " both empty", 1.0, empty, EPS);
        }
    }

    /** Containment is never below Jaccard, since the smaller set is never larger than the union. */
    @Test
    public void testContainmentIsNeverBelowJaccard() {
        long[][] cases = { { 10, 100, 20 }, { 0, 5, 900 }, { 7, 7, 7 }, { 1, 1_000_000, 2 } };
        for (long[] c : cases) {
            assertTrue("shared=" + c[0],
                    similarity(CONTAINMENT, c[0], c[1], c[2])
                            >= similarity(JACCARD, c[0], c[1], c[2]) - EPS);
        }
    }

    /** Jaccard is what a project gets without saying anything; {@code false} is the way to containment. */
    @Test
    public void testTheConfigKeyDefaultsToJaccard() {
        ConfigParamInfo<?> info = FTConfigKey.JACCARD_SIM.getInfo();
        assertEquals("jaccardSim", FTConfigKey.JACCARD_SIM.getName());
        assertEquals(Boolean.TRUE, info.defaultValue());
        assertTrue(info.isValid("true"));
        assertTrue(info.isValid("false"));
        assertTrue(info.isValid(info.getMDDefaultValue()));
    }
}
