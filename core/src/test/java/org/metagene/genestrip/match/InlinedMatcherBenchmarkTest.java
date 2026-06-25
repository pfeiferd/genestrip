/*
 *
 * "Commons Clause" License Condition v1.0
 *
 * The Software is provided to you by the Licensor under the License,
 * as defined below, subject to the following condition.
 *
 * Without limiting other conditions in the License, the grant of rights under the License
 * will not include, and the License does not grant to you, the right to Sell the Software.
 *
 * For purposes of the foregoing, "Sell" means practicing any or all of the rights granted
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
package org.metagene.genestrip.match;

import it.unimi.dsi.fastutil.objects.Object2LongMap;
import org.junit.Test;
import org.metagene.genestrip.DefaultExecutionContext;
import org.metagene.genestrip.ExecutionContext;
import org.metagene.genestrip.match.FastqKMerMatcher.MatcherReadEntry;
import org.metagene.genestrip.store.KMerSortedArray;
import org.metagene.genestrip.store.KMerSortedArray.ValueConverter;
import org.metagene.genestrip.store.KMerUniqueCounterBits;
import org.metagene.genestrip.tax.SmallTaxTree.SmallTaxIdNode;
import org.metagene.genestrip.util.CGAT;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

import static org.junit.Assert.*;

/**
 * Regression test: verifies that FastqKMerMatcher and InlinedFastqKMerMatcher
 * produce byte-for-byte identical statistics on identical input.
 *
 * Benchmark test: measures the speedup factor of the inlined variant.
 * Uses k=31, realistic read length, and a sorted k-mer store to reflect
 * production workloads. Runs both matchers with a JIT warm-up phase before
 * timing, then reports the speedup and an estimated performance range.
 */
public class InlinedMatcherBenchmarkTest {

    // --- Configuration -------------------------------------------------------

    // The regression and benchmark are repeated for each of these k-mer sizes.
    protected static final int[] K_VALUES = {16, 21, 31};
    // Three taxonomy ids spread across the store.
    protected static final String[] TAXIDS = {"1", "2", "3"};
    // k-mers per taxid for the performance benchmark (store sized for a realistic workload).
    protected static final int KMERS_PER_TAXID = 50000000;
    // k-mers per taxid for the correctness regression — kept small (1500 total) so the
    // matcher-equivalence check stays fast and cheap regardless of the benchmark store size.
    protected static final int REGRESSION_KMERS_PER_TAXID = 500;
    // Typical Illumina short-read length.
    protected static final int READ_LENGTH = 250;
    // How many known k-mers to embed per read (controls hit density).
    protected static final int KMERS_EMBEDDED_PER_READ = (int) (READ_LENGTH * 0.2);

    // --- Setup helpers -------------------------------------------------------

    /**
     * Builds a sorted store of the given k-mer size {@code k} with {@code kmersPerTaxid}
     * randomly generated k-mers per taxid. Stores the straight-encoded bytes of each
     * canonical k-mer in the returned setup so they can be embedded in test reads later.
     */
    @SuppressWarnings("unchecked")
    static TestSetup buildSetup(int k, int kmersPerTaxid, Random rng) {
        int totalKmers = kmersPerTaxid * TAXIDS.length;
        KMerSortedArray<String> strStore = new KMerSortedArray<>(
                k, 0.001, 0.001, Arrays.asList(TAXIDS), false, true);
        strStore.initSize(totalKmers);

        List<byte[]> allKmers = new ArrayList<>(totalKmers);
        byte[] seq = new byte[k];
        for (String taxid : TAXIDS) {
            int added = 0;
            while (added < kmersPerTaxid) {
                for (int j = 0; j < k; j++) {
                    seq[j] = CGAT.DECODE_TABLE[rng.nextInt(4)];
                }
                long canonical = CGAT.kMerToLong(seq, 0, k, null);
                if (canonical == -1) {
                    continue; // bad char from kMerToLong (should not happen with DECODE_TABLE)
                }
                if (strStore.putLong(canonical, taxid)) {
                    // Store the byte form of the canonical k-mer for embedding in reads.
                    // longToKMerStraight is the exact inverse of kMerToLongStraight, so
                    // when the matcher re-encodes these bytes it will recover `canonical`.
                    byte[] kmerBytes = new byte[k];
                    CGAT.longToKMerStraight(canonical, kmerBytes, 0, k);
                    allKmers.add(kmerBytes);
                    added++;
                }
            }
        }
        strStore.optimize();

        // Convert String -> SmallTaxIdNode, preserving storeIndex assignments.
        KMerSortedArray<SmallTaxIdNode> nodeStore = new KMerSortedArray<>(strStore,
                new ValueConverter<String, SmallTaxIdNode>() {
                    @Override
                    public SmallTaxIdNode convertValue(String value) {
                        SmallTaxIdNode node = new SmallTaxIdNode(value, null, null);
                        node.setStoreIndex(strStore.getIndexForValue(value));
                        return node;
                    }
                });

        return new TestSetup(strStore, nodeStore, allKmers);
    }

    /**
     * Generates {@code count} reads of length READ_LENGTH.
     * Each read is random DNA with KMERS_EMBEDDED_PER_READ known k-mers
     * from {@code allKmers} embedded at random positions.
     */
    static byte[][] generateReads(int k, Random rng, List<byte[]> allKmers, int count) {
        byte[][] reads = new byte[count][READ_LENGTH];
        for (int i = 0; i < count; i++) {
            byte[] read = reads[i];
            for (int j = 0; j < READ_LENGTH; j++) {
                read[j] = CGAT.DECODE_TABLE[rng.nextInt(4)];
            }
            for (int hit = 0; hit < KMERS_EMBEDDED_PER_READ; hit++) {
                int pos = rng.nextInt(READ_LENGTH - k + 1);
                byte[] kmer = allKmers.get(rng.nextInt(allKmers.size()));
                System.arraycopy(kmer, 0, read, pos, k);
            }
        }
        return reads;
    }

    /**
     * Creates a single-threaded matcher (taxTree=null, no classification path).
     */
    static FastqKMerMatcher createMatcher(boolean inlined, KMerSortedArray<SmallTaxIdNode> store) {
        ExecutionContext bundle = new DefaultExecutionContext(null, 0, 1000);
        if (inlined) {
            return new InlinedFastqKMerMatcher(store, READ_LENGTH * 2, 1000,
                    bundle, false, 0, null, 4, -1, -1, true, 1, null);
        }
        return new FastqKMerMatcher(store, READ_LENGTH * 2, 1000,
                bundle, false, 0, null, 4, -1, -1, true, 1, null);
    }

    /**
     * Runs all reads through the matcher, resetting the entry completely
     * for each read (regression-safe path).
     */
    static void runReadsForStats(FastqKMerMatcher matcher, MatcherReadEntry entry,
            byte[][] reads, long startReadNo) {
        entry.readNo = startReadNo;
        for (byte[] read : reads) {
            System.arraycopy(read, 0, entry.read, 0, READ_LENGTH);
            entry.readSize = READ_LENGTH;
            entry.readNo++;
            entry.bufferPos = 0;
            entry.usedPaths = 0;
            entry.classNode = null;
            Arrays.fill(entry.readTaxIdNode, null);
            Arrays.fill(entry.counts, (short) 0);
            matcher.matchRead(entry, 0);
        }
    }

    /**
     * Runs all reads through the matcher with minimal overhead reset.
     * With taxTree=null, classification state (classNode, usedPaths, …) is never
     * written by matchRead, so only readNo, readSize, and read data matter.
     */
    static void runReadsForBenchmark(FastqKMerMatcher matcher, MatcherReadEntry entry,
            byte[][] reads) {
        for (byte[] read : reads) {
            System.arraycopy(read, 0, entry.read, 0, READ_LENGTH);
            entry.readSize = READ_LENGTH;
            entry.readNo++;
            matcher.matchRead(entry, 0);
        }
    }

    // --- Regression test -----------------------------------------------------

    /**
     * Verifies that every statistics field produced by FastqKMerMatcher and
     * InlinedFastqKMerMatcher is identical when given the same reads in the
     * same order. Repeated for each k-mer size in {@link #K_VALUES}.
     *
     * Fields checked per taxid:
     *   kmers, reads1KMer, contigs, maxContigLen, contigLenSquaredSum
     * Plus: unique k-mer counts via KMerUniqueCounterBits.
     */
    @Test
    public void testIdenticalResults() {
        for (int k : K_VALUES) {
            TestSetup setup = buildSetup(k, REGRESSION_KMERS_PER_TAXID, new Random(12345L));
            byte[][] reads = generateReads(k, new Random(67890L), setup.allKmers, 500);
            KMerUniqueCounterBits uniqueStdCounter = new KMerUniqueCounterBits(setup.strStore, false);
            KMerUniqueCounterBits uniqueInlCounter = new KMerUniqueCounterBits(setup.strStore, false);

            // --- Standard matcher ---
            FastqKMerMatcher stdMatcher = createMatcher(false, setup.nodeStore);
            MatcherReadEntry stdEntry = new MatcherReadEntry(READ_LENGTH * 2, false, 4);
            stdMatcher.initStats();
            stdMatcher.initUniqueCounter(uniqueStdCounter);
            runReadsForStats(stdMatcher, stdEntry, reads, 0L);
            Object2LongMap<String> uniqueStd = uniqueStdCounter.getUniqueKmerCounts();

            // --- Inlined matcher on identical reads ---
            InlinedFastqKMerMatcher inlMatcher =
                    (InlinedFastqKMerMatcher) createMatcher(true, setup.nodeStore);
            MatcherReadEntry inlEntry = new MatcherReadEntry(READ_LENGTH * 2, false, 4);
            inlMatcher.initStats();
            inlMatcher.initUniqueCounter(uniqueInlCounter);
            runReadsForStats(inlMatcher, inlEntry, reads, 0L);
            Object2LongMap<String> uniqueInl = uniqueInlCounter.getUniqueKmerCounts();

            // --- Compare all statistics per taxid ---
            for (String taxid : TAXIDS) {
                int vi = setup.strStore.getIndexForValue(taxid);
                CountsPerTaxid statsStd = stdMatcher.statsIndex[vi];
                CountsPerTaxid statsInl = inlMatcher.statsIndex[vi];

                if (statsStd == null) {
                    assertNull("k=" + k + " taxid=" + taxid + ": inlined unexpectedly non-null", statsInl);
                    continue;
                }
                assertNotNull("k=" + k + " taxid=" + taxid + ": inlined unexpectedly null", statsInl);

                assertEquals("k=" + k + " taxid=" + taxid + " kmers",
                        statsStd.getKMers(), statsInl.getKMers());
                assertEquals("k=" + k + " taxid=" + taxid + " reads1KMer",
                        statsStd.getReads1KMer(), statsInl.getReads1KMer());
                assertEquals("k=" + k + " taxid=" + taxid + " contigs",
                        statsStd.getContigs(), statsInl.getContigs());
                assertEquals("k=" + k + " taxid=" + taxid + " maxContigLen",
                        statsStd.getMaxContigLen(), statsInl.getMaxContigLen());
                assertEquals("k=" + k + " taxid=" + taxid + " contigLenSquaredSum",
                        statsStd.contigLenSquaredSum, statsInl.contigLenSquaredSum);

                // Unique k-mer counts (position-based via KMerUniqueCounterBits)
                assertEquals("k=" + k + " taxid=" + taxid + " uniqueKmers",
                        uniqueStd.getLong(taxid), uniqueInl.getLong(taxid));
            }
        }
    }

    // --- Performance benchmark -----------------------------------------------

    /**
     * Benchmarks FastqKMerMatcher vs InlinedFastqKMerMatcher with realistic read length,
     * repeated for each k-mer size in {@link #K_VALUES}.
     *
     * Protocol (per k):
     *  1. JIT warm-up: 20 000 reads through each matcher before timing starts.
     *  2. 7 timed rounds of 10 000 reads each, run alternately per matcher
     *     to spread any OS-level jitter equally.
     *  3. The fastest and slowest round are discarded; the remaining 5 are averaged.
     *
     * Prints timing and speedup factor to stdout.
     * Asserts speedup >= 0.7 as a sanity bound
     * (the inlined version should not be dramatically slower).
     *
     * Realistic expected range: 1.05x – 1.25x speedup, primarily from:
     *  - Caching CGAT lookup tables as local finals (avoids static field loads)
     *  - Reading entry.read[i+k-1] only once per increment step
     *  - Enabling the JIT to see the entire hot loop without indirect calls
     */
    @Test
    public void benchmarkSpeedup() {
        final int WARMUP_READS = 20_000;
        final int BENCH_ROUNDS = 7;       // discard min + max → 5 effective rounds
        final int READS_PER_ROUND = 10_000;

        for (int k : K_VALUES) {
            TestSetup setup = buildSetup(k, KMERS_PER_TAXID, new Random(99999L));
            int totalReads = WARMUP_READS + READS_PER_ROUND * BENCH_ROUNDS;
            byte[][] reads = generateReads(k, new Random(11111L), setup.allKmers, totalReads);

            FastqKMerMatcher stdMatcher = createMatcher(false, setup.nodeStore);
            InlinedFastqKMerMatcher inlMatcher =
                    (InlinedFastqKMerMatcher) createMatcher(true, setup.nodeStore);
            MatcherReadEntry stdEntry = new MatcherReadEntry(READ_LENGTH * 2, false, 4);
            MatcherReadEntry inlEntry = new MatcherReadEntry(READ_LENGTH * 2, false, 4);
            stdEntry.readNo = 0;
            inlEntry.readNo = 0;

            // --- Warm-up ---
            stdMatcher.initStats();
            inlMatcher.initStats();
            byte[][] warmupReads = Arrays.copyOfRange(reads, 0, WARMUP_READS);
            runReadsForBenchmark(stdMatcher, stdEntry, warmupReads);
            runReadsForBenchmark(inlMatcher, inlEntry, warmupReads);
            // Reset stats to avoid skew from accumulated entries during warm-up.
            stdMatcher.initStats();
            inlMatcher.initStats();

            // --- Timed rounds (alternating to share OS jitter fairly) ---
            long[] stdNs = new long[BENCH_ROUNDS];
            long[] inlNs = new long[BENCH_ROUNDS];
            for (int r = 0; r < BENCH_ROUNDS; r++) {
                byte[][] roundReads = Arrays.copyOfRange(reads,
                        WARMUP_READS + r * READS_PER_ROUND,
                        WARMUP_READS + (r + 1) * READS_PER_ROUND);

                long t0 = System.nanoTime();
                runReadsForBenchmark(stdMatcher, stdEntry, roundReads);
                stdNs[r] = System.nanoTime() - t0;

                long t1 = System.nanoTime();
                runReadsForBenchmark(inlMatcher, inlEntry, roundReads);
                inlNs[r] = System.nanoTime() - t1;
            }

            // Discard min + max, average the 5 middle rounds.
            Arrays.sort(stdNs);
            Arrays.sort(inlNs);
            long stdSum = 0, inlSum = 0;
            int effective = BENCH_ROUNDS - 2;
            for (int r = 1; r <= effective; r++) {
                stdSum += stdNs[r];
                inlSum += inlNs[r];
            }
            double stdMs = stdSum / (1_000_000.0 * effective);
            double inlMs = inlSum / (1_000_000.0 * effective);
            double speedup = stdMs / inlMs;

            long innerLoopsPerRound = (long) READS_PER_ROUND * (READ_LENGTH - k + 1);
            System.out.printf("%n=== InlinedFastqKMerMatcher Benchmark (k=%d) ===%n", k);
            System.out.printf("  Store:  k=%d, %d k-mers (%d taxids x %d)%n",
                    k, KMERS_PER_TAXID * TAXIDS.length, TAXIDS.length, KMERS_PER_TAXID);
            System.out.printf("  Reads:  %d reads x %d bp per round (%,d inner-loop iterations)%n",
                    READS_PER_ROUND, READ_LENGTH, innerLoopsPerRound);
            System.out.printf("  Rounds: %d timed, min+max discarded, %d averaged%n",
                    BENCH_ROUNDS, effective);
            System.out.printf("  FastqKMerMatcher (standard):  %6.2f ms/round  (%5.1f ns/read)%n",
                    stdMs, stdMs * 1_000_000 / READS_PER_ROUND);
            System.out.printf("  InlinedFastqKMerMatcher:       %6.2f ms/round  (%5.1f ns/read)%n",
                    inlMs, inlMs * 1_000_000 / READS_PER_ROUND);
            System.out.printf("  Speedup factor: %.3fx%n", speedup);
            if (speedup >= 1.0) {
                System.out.printf("  -> Inlined is %.1f%% faster%n", (speedup - 1.0) * 100);
            } else {
                System.out.printf("  -> Inlined is %.1f%% slower (JIT may already inline standard variant)%n",
                        (1.0 - speedup) * 100);
            }

            // Sanity: inlined should not be dramatically slower
            assertTrue("Speedup should be >= 0.7 for k=" + k + " (got " + String.format("%.3f", speedup) + ")",
                    speedup >= 0.7);
        }
    }

    // --- Shared data holder --------------------------------------------------

    static class TestSetup {
        final KMerSortedArray<String> strStore;
        final KMerSortedArray<SmallTaxIdNode> nodeStore;
        final List<byte[]> allKmers;

        TestSetup(KMerSortedArray<String> strStore,
                KMerSortedArray<SmallTaxIdNode> nodeStore,
                List<byte[]> allKmers) {
            this.strStore = strStore;
            this.nodeStore = nodeStore;
            this.allKmers = allKmers;
        }
    }
}
