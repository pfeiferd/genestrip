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
package org.metagene.genestrip.finertree.probfilter;

import org.junit.Test;
import org.metagene.genestrip.probfilter.BlockedBloomFilter;
import org.metagene.genestrip.probfilter.BloomFilter;
import org.metagene.genestrip.probfilter.MurmurBloomFilter;
import org.metagene.genestrip.probfilter.ProbFilter;
import org.metagene.genestrip.probfilter.SingleWordBloomFilter;
import org.metagene.genestrip.probfilter.XORBlockedBloomFilter;
import org.metagene.genestrip.probfilter.XORBloomFilter;
import org.metagene.genestrip.probfilter.XORSingleWordBloomFilter;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

/**
 * Tests the properties of {@link KMerIndexFilterHelper#combine(long, int)} that the filters it feeds
 * rely on, and checks the outcome against every filter of the {@code probfilter} package.
 * <p>
 * A filter takes a key's bit positions from the key folded by a 32 bit rotation, and its word position
 * from the hash. An index that ends up in a rotation invariant mask cancels in that fold, and one that
 * stays in the low bits is dropped again wherever the word position comes from the upper half of a
 * product. A k-mer stored under one index would then be reported as present under every other one -
 * which is what {@code KMerStoreWorkGoal} asks the filter about for every child of a refinement node.
 * <p>
 * No filter of the package is currently exposed to that, which is worth stating so the sweep in
 * {@link #testNoCrossIndexFalsePositives()} is not read as more than it is. The mixing filters put a
 * MurmurHash3 finalizer in front of the fold, so every distinct key reaches them; the {@code XOR}
 * variants do hash trivially as {@code seed ^ key}, but they reduce by a modulo over every bit of it -
 * the pairing {@code HashReducePairingTest} pins. Only a non-mixing hash left on a multiply-shift
 * reduction is vulnerable, and the package no longer contains one. Measured against such a filter a
 * rotation invariant index yields 86% cross index false positives, against 1.3% for {@code combine};
 * against every filter as it stands today the two are indistinguishable. The sweep is therefore a
 * regression guard for a filter that reintroduces the hazard, and the structural tests below are what
 * actually pin {@code combine}.
 */
public class KMerIndexFilterHelperTest {
    /**
     * Number of store indexes covered. They are bounded by the store's {@code MAX_VALUES}; a viral
     * database has roughly 46000 tree nodes and hence that many distinct ones.
     */
    private static final int INDEXES = 65535;

    /** An arbitrary but fixed k-mer, so that all tests here are deterministic. */
    private static final long KMER = new Random(7).nextLong() >>> 2;

    /** Keys entered into each filter of {@link #testNoCrossIndexFalsePositives()}. */
    private static final int ENTRIES = 200000;
    /** Lookups under a foreign index performed per entered k-mer. */
    private static final int PROBES_PER_KMER = 10;
    /** Distinct store indexes drawn from, i.e. roughly the tree node count of a viral database. */
    private static final int INDEX_RANGE = 46000;
    /**
     * Memory each filter is given, so that the rates below are comparable. The classical family sizes
     * itself from a false-positive probability instead, for which {@link #FPP} is the closest match.
     */
    private static final int BITS_PER_KEY = 10;
    /** Target false-positive probability of the classical family, i.e. about {@link #BITS_PER_KEY}. */
    private static final double FPP = 0.01;
    /**
     * Rate every filter has to stay below. Well above what any of them produces on its own at
     * {@link #BITS_PER_KEY} bits per key (under 2%), and well below the 86% a filter whose word
     * position ignores the index produces.
     */
    private static final double MAX_CROSS_INDEX_RATE = 0.05;

    /**
     * Combining one k-mer with different indexes must yield different keys - otherwise the pair cannot
     * be told apart at all. {@code Integer.MAX_VALUE} is included because it is the marker
     * {@code KMerIndexBloomGoal.OTHER_VALUE} uses for the "OTHER" bucket.
     */
    @Test
    public void testDistinctKeyPerIndex() {
        Set<Long> keys = new HashSet<>();
        for (int i = 0; i < INDEXES; i++) {
            keys.add(KMerIndexFilterHelper.combine(KMER, i));
        }
        assertEquals(INDEXES, keys.size());
        assertTrue(keys.add(KMerIndexFilterHelper.combine(KMER, Integer.MAX_VALUE)));
    }

    /**
     * The index must survive the 32 bit rotation fold the filters derive a key's bit positions from.
     * Folding cancels every rotation invariant part of a key, so an index xored once into each half of
     * the word would collapse all indexes onto a single set of bit positions.
     */
    @Test
    public void testIndexSurvivesRotationFold() {
        Set<Long> folded = new HashSet<>();
        for (int i = 0; i < INDEXES; i++) {
            long key = KMerIndexFilterHelper.combine(KMER, i);
            folded.add(key ^ Long.rotateLeft(key, 32));
        }
        // The fold maps onto 2^32 values, so a few birthday collisions among 65535 indexes are expected
        // and fine. Collapsing onto one value is not.
        assertTrue("distinct folded keys: " + folded.size(), folded.size() > INDEXES * 0.99);
    }

    /**
     * The index must also reach the upper bits, because the word position is taken from the upper half
     * of a product. An index xored in unmixed only occupies low bits, which are dropped there.
     */
    @Test
    public void testIndexReachesHighBits() {
        Set<Long> high = new HashSet<>();
        for (int i = 0; i < INDEXES; i++) {
            high.add(KMerIndexFilterHelper.combine(KMER, i) >>> 48);
        }
        // 65535 indexes drawn over 16 bits fill about 1 - 1/e of them; anything near 1 means the index
        // does not reach the upper bits at all.
        assertTrue("distinct high bits: " + high.size(), high.size() > 30000);
    }

    /**
     * End to end check against the filters themselves: every k-mer is entered under exactly one index
     * and then looked up under other ones, all of which must report absence apart from the filter's own
     * false positives.
     * <p>
     * Run against every filter of the package, since {@code combine} has to serve whichever one
     * {@code KMerIndexBloomGoal} is built with. As the class comment records, none of them can currently
     * tell {@code combine} apart from a lazier folding of the index, so this passing says the pair is
     * sound for the filters that exist, not that {@code combine} is what makes it so. Every filter is
     * measured before anything is asserted, so that one bad variant does not hide the others.
     */
    @Test
    public void testNoCrossIndexFalsePositives() {
        StringBuilder measured = new StringBuilder();
        List<String> tooHigh = new ArrayList<>();
        for (Map.Entry<String, ProbFilter> entry : newFilters().entrySet()) {
            double rate = crossIndexFalsePositiveRate(entry.getValue());
            measured.append(' ').append(entry.getKey()).append('=').append(rate);
            if (rate >= MAX_CROSS_INDEX_RATE) {
                tooHigh.add(entry.getKey());
            }
        }
        assertTrue("cross index false positive rate at or above " + MAX_CROSS_INDEX_RATE + " for " + tooHigh
                + " - measured:" + measured, tooHigh.isEmpty());
    }

    /**
     * Returns every concrete {@link ProbFilter} of the package, each given {@link #BITS_PER_KEY} bits
     * per key so that the measured rates are comparable. A new implementation belongs here, or
     * {@code combine} is never checked against the way it derives its bit and word positions.
     *
     * @return the filters to measure, keyed by the name reported on failure
     */
    private Map<String, ProbFilter> newFilters() {
        Map<String, ProbFilter> filters = new LinkedHashMap<>();
        filters.put("BloomFilter", new BloomFilter(FPP, ENTRIES));
        filters.put("MurmurBloomFilter", new MurmurBloomFilter(FPP, ENTRIES));
        filters.put("XORBloomFilter", new XORBloomFilter(FPP, ENTRIES));
        filters.put("BlockedBloomFilter", new BlockedBloomFilter(ENTRIES, BITS_PER_KEY));
        filters.put("XORBlockedBloomFilter", new XORBlockedBloomFilter(ENTRIES, BITS_PER_KEY));
        filters.put("SingleWordBloomFilter", new SingleWordBloomFilter(ENTRIES, BITS_PER_KEY));
        filters.put("XORSingleWordBloomFilter", new XORSingleWordBloomFilter(ENTRIES, BITS_PER_KEY));
        return filters;
    }

    /**
     * Fills the given filter with {@link #ENTRIES} (k-mer, index) pairs and returns the fraction of
     * lookups under a <em>different</em> index that it nevertheless reports as present. The k-mers and
     * indexes are drawn from fixed seeds, so every filter sees the same keys and the result is
     * deterministic.
     *
     * @param filter the filter to fill and measure
     * @return the measured cross index false-positive rate
     */
    private double crossIndexFalsePositiveRate(ProbFilter filter) {
        Random random = new Random(999);
        long[] kmers = new long[ENTRIES];
        int[] storeIndexes = new int[ENTRIES];
        for (int i = 0; i < ENTRIES; i++) {
            kmers[i] = random.nextLong() >>> 2;
            storeIndexes[i] = random.nextInt(INDEX_RANGE);
            filter.putLong(KMerIndexFilterHelper.combine(kmers[i], storeIndexes[i]));
        }

        long falsePositives = 0;
        long lookups = 0;
        Random probes = new Random(4242);
        for (int i = 0; i < ENTRIES; i++) {
            for (int p = 0; p < PROBES_PER_KMER; p++) {
                int other = probes.nextInt(INDEX_RANGE);
                if (other == storeIndexes[i]) {
                    continue;
                }
                lookups++;
                if (filter.containsLong(KMerIndexFilterHelper.combine(kmers[i], other))) {
                    falsePositives++;
                }
            }
        }
        return (double) falsePositives / lookups;
    }
}
