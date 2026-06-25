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

import org.junit.Test;
import org.metagene.genestrip.store.KMerSortedArray;
import org.metagene.genestrip.store.KMerSortedArray.ValueConverter;
import org.metagene.genestrip.tax.SmallTaxTree.SmallTaxIdNode;
import org.metagene.genestrip.util.CGAT;

import java.util.Arrays;
import java.util.Random;

import static org.junit.Assert.*;

/**
 * Correctness test: verifies that the radix-guided search and binary search return
 * identical results for any query — both for stored k-mers (hits) and for
 * random k-mers not in the store (misses).
 *
 * Benchmark test: measures binary search vs radix search by directly
 * calling {@code getLong} in a tight loop on a store large enough to exceed
 * L3 cache. Uses {@code enforceLarge=true} to exercise the BigArray code path
 * that matters for production databases.
 *
 * <p>The radix search (see {@link KMerSortedArray#setRadixSearch(boolean)}) takes the
 * next 4 key bits at each step and probes into the corresponding one of 16 equal
 * sub-ranges of the current interval (a single {@code >>> 4} shift), shrinking the
 * interval by ~16x per step instead of 2x. The expectation — the reason it may beat
 * both binary and interpolation search out of cache — is that the early probes land on
 * a small fixed set of array positions that stay warm in the L3 cache, while the total
 * probe count drops from ~log2(n) to ~log2(n)/4.
 *
 * No lists of stored k-mers are kept — the sorted store is sampled in O(1)
 * per query via {@link KMerSortedArray#getKMerAt(long)}.
 */
public class RadixSearchBenchmarkTest extends InlinedMatcherBenchmarkTest {

    // K and TAXIDS are reused from InlinedMatcherBenchmarkTest (k=31, three taxids); only
    // the search is exercised here, so the exact taxid set is irrelevant beyond being non-empty.

    // --- Configuration -------------------------------------------------------

    // Radix widths (bits consumed per step → 2^bits buckets) swept by the benchmark and
    // verified for correctness. 4 is the default; larger widths do fewer but pricier steps.
    private static final int[] RADIX_BITS_VALUES = {4, 6, 8, 10, 12};

    // Medium store for correctness test (fast to build, all entries verified).
    private static final int CORRECTNESS_KMERS_PER_TAXID = 20_000; // 60 K total

    // Large store for benchmark. enforceLarge=true forces the BigArray path.
    // At 8 bytes/entry the default 500 M × 3 = 1.5 B entries → ~12 GB of k-mer data, which
    // is the out-of-cache regime where the reduced probe count is expected to pay off.
    // Override with -Dgenestrip.bench.kmersPerTaxid=<n> to run the benchmark on a smaller
    // heap (e.g. a few million) without recompiling.
    private static final int BENCH_KMERS_PER_TAXID =
            Integer.getInteger("genestrip.bench.kmersPerTaxid", 500_000_00);

    // Number of k-mer queries used per benchmark round (hits + misses, 50/50).
    // Only this many k-mers are held in memory at a time — independent of store size.
    private static final int QUERIES_PER_ROUND = 100_000;

    // --- Store builder -------------------------------------------------------

    /**
     * Builds a sorted store of the given k-mer size {@code k} from randomly generated k-mers.
     * No list of stored k-mers is maintained; use {@link #sampleKMers} to
     * obtain representative query k-mers after building.
     *
     * @param k             k-mer size
     * @param kmersPerTaxid number of k-mers to insert per taxid
     * @param enforceLarge  when {@code true}, forces use of the BigArray storage
     *                      path regardless of the entry count
     * @param rng           random source for k-mer generation
     * @return converted SmallTaxIdNode store, sorted and ready for queries
     */
    private static KMerSortedArray<SmallTaxIdNode> buildStore(
            int k, int kmersPerTaxid, boolean enforceLarge, Random rng) {

        long total = (long) kmersPerTaxid * TAXIDS.length;
        KMerSortedArray<String> strStore = new KMerSortedArray<>(
                k, 0.001, 0.001, Arrays.asList(TAXIDS), enforceLarge, true);
        strStore.initSize(total);

        byte[] seq = new byte[k];
        for (String taxid : TAXIDS) {
            int added = 0;
            while (added < kmersPerTaxid) {
                for (int j = 0; j < k; j++) seq[j] = CGAT.DECODE_TABLE[rng.nextInt(4)];
                long canonical = CGAT.kMerToLong(seq, 0, k, null);
                if (canonical != -1 && strStore.putLong(canonical, taxid)) {
                    added++;
                }
            }
        }
        strStore.optimize();

        return new KMerSortedArray<>(strStore, new ValueConverter<String, SmallTaxIdNode>() {
            @Override
            public SmallTaxIdNode convertValue(String value) {
                SmallTaxIdNode node = new SmallTaxIdNode(value, null, null);
                node.setStoreIndex(strStore.getIndexForValue(value));
                return node;
            }
        });
    }

    /**
     * Samples {@code n} k-mers from uniformly random positions in the sorted
     * store. O(n) time, O(n) space — independent of store size.
     */
    private static long[] sampleKMers(KMerSortedArray<?> store, int n, Random rng) {
        long entries = store.getEntries();
        long[] sample = new long[n];
        for (int i = 0; i < n; i++) {
            // nextLong() >>> 1 yields a non-negative long; % entries maps it to [0, entries).
            sample[i] = store.getKMerAt((rng.nextLong() >>> 1) % entries);
        }
        return sample;
    }

    /**
     * Generates {@code n} random k-mers that are (with overwhelming probability)
     * not in any store.  Used as miss queries.
     */
    private static long[] generateMissKMers(int k, int n, Random rng) {
        long[] misses = new long[n];
        byte[] seq = new byte[k];
        for (int i = 0; i < n; i++) {
            long canonical;
            do {
                for (int j = 0; j < k; j++) seq[j] = CGAT.DECODE_TABLE[rng.nextInt(4)];
                canonical = CGAT.kMerToLong(seq, 0, k, null);
            } while (canonical == -1);
            misses[i] = canonical;
        }
        return misses;
    }

    // --- Correctness test ----------------------------------------------------

    /** Regular {@code long[]} array path. */
    @Test
    public void testRadixIdenticalToBinarySmall() {
        assertRadixIdenticalToBinary(false);
    }

    /** Large {@code long[][]} BigArray path (different probe index arithmetic). */
    @Test
    public void testRadixIdenticalToBinaryLarge() {
        assertRadixIdenticalToBinary(true);
    }

    /**
     * For a sample of stored k-mers and a large set of random (mostly non-stored)
     * k-mers, asserts that binary search and radix search return identical
     * results and identical posStore values for every query.
     */
    private void assertRadixIdenticalToBinary(boolean enforceLarge) {
        for (int k : K_VALUES) {
            KMerSortedArray<SmallTaxIdNode> store =
                    buildStore(k, CORRECTNESS_KMERS_PER_TAXID, enforceLarge, new Random(54321L));
            // Disable the Bloom filter so misses actually reach the search and are compared.
            store.setUseFilter(false);

            long[] posStoreBinary = new long[1];
            long[] posStoreRadix = new long[1];
            long totalEntries = store.getEntries();

            // Verify every benchmarked radix width against binary search.
            for (int bits : RADIX_BITS_VALUES) {
                store.setRadixBits(bits);
                Random rng = new Random(98765L); // identical miss queries across widths

                // --- Hits: sample all stored k-mers from the sorted array ---
                int hitCount = 0;
                for (long pos = 0; pos < totalEntries; pos++) {
                    long kmer = store.getKMerAt(pos);

                    store.setRadixSearch(false);
                    SmallTaxIdNode binaryResult = store.getLong(kmer, posStoreBinary);

                    store.setRadixSearch(true);
                    SmallTaxIdNode radixResult = store.getLong(kmer, posStoreRadix);

                    assertNotNull("k=" + k + " bits=" + bits + " binary search missed stored kmer at pos " + pos, binaryResult);
                    assertNotNull("k=" + k + " bits=" + bits + " radix search missed stored kmer at pos " + pos, radixResult);
                    assertSame("k=" + k + " bits=" + bits + " different node for stored kmer at pos " + pos,
                            binaryResult, radixResult);
                    assertEquals("k=" + k + " bits=" + bits + " posStore differs for stored kmer at pos " + pos,
                            posStoreBinary[0], posStoreRadix[0]);
                    hitCount++;
                }

                // --- Random queries (mostly misses) ---
                int missCount = 0, randomHits = 0;
                byte[] seq = new byte[k];
                for (int i = 0; i < 100_000; i++) {
                    for (int j = 0; j < k; j++) seq[j] = CGAT.DECODE_TABLE[rng.nextInt(4)];
                    long kmer = CGAT.kMerToLong(seq, 0, k, null);
                    if (kmer == -1) continue;

                    store.setRadixSearch(false);
                    SmallTaxIdNode binaryResult = store.getLong(kmer, posStoreBinary);

                    store.setRadixSearch(true);
                    SmallTaxIdNode radixResult = store.getLong(kmer, posStoreRadix);

                    assertEquals("k=" + k + " bits=" + bits + " binary and radix disagree on random kmer " + i,
                            binaryResult == null, radixResult == null);
                    if (binaryResult != null) {
                        assertSame("k=" + k + " bits=" + bits + " different node for random kmer " + i,
                                binaryResult, radixResult);
                        assertEquals("k=" + k + " bits=" + bits + " posStore differs for random kmer " + i,
                                posStoreBinary[0], posStoreRadix[0]);
                        randomHits++;
                    } else {
                        missCount++;
                    }
                }

                System.out.printf("%nCorrectness (%s path, k=%d, radix=%d bits): %d hits, %d misses, %d random hits verified%n",
                        enforceLarge ? "large/BigArray" : "small", k, bits, hitCount, missCount, randomHits);
            }
        }
    }

    // --- Benchmark -----------------------------------------------------------

    /**
     * Benchmarks binary search against the radix search at several radix widths, repeated for
     * each k-mer size in {@link #K_VALUES}. For every width in {@link #RADIX_BITS_VALUES} it
     * times both the plain lookup ({@link KMerSortedArray#getLong}) and the combined variant
     * with the manually inlined lookup ({@link KMerSortedArray#getLongInlined}), all relative
     * to the binary-search baseline (radix off).
     *
     * <p>The store is built with {@code enforceLarge=true} to exercise the
     * BigArray code path. Only {@value #QUERIES_PER_ROUND} query k-mers
     * are held in memory at any time (sampled from the store), so memory
     * usage is independent of store size beyond the store itself.
     *
     * <p>Query mix: 50% hits (k-mers sampled from the store), 50% misses
     * (random k-mers). The Bloom filter short-circuits most misses; the
     * binary/radix search is exercised mainly for hits.
     */
    @Test
    public void benchmarkRadixVsBinary() {
        final int WARMUP_ITERS = 5;
        final int BENCH_ROUNDS = 7;
        final int HALF = QUERIES_PER_ROUND / 2;
        final int nWidths = RADIX_BITS_VALUES.length;

        for (int k : K_VALUES) {
            Random rng = new Random(77777L);

            // Build with enforceLarge=true to test the BigArray (large store) code path.
            KMerSortedArray<SmallTaxIdNode> store =
                    buildStore(k, BENCH_KMERS_PER_TAXID, true, rng);

            long entries = store.getEntries();

            // Sample hit and miss queries — small fixed arrays, independent of store size.
            long[] hits   = sampleKMers(store, HALF, rng);
            long[] misses = generateMissKMers(k, HALF, rng);

            // Interleave hits and misses for a realistic query stream.
            long[] queries = new long[QUERIES_PER_ROUND];
            for (int i = 0; i < HALF; i++) {
                queries[2 * i]     = hits[i];
                queries[2 * i + 1] = misses[i];
            }

            long[] posStore = new long[1];

            // --- Warm-up (all variants and widths, several passes) ---
            for (int w = 0; w < WARMUP_ITERS; w++) {
                store.setRadixSearch(false);
                for (long q : queries) store.getLong(q, posStore);

                store.setRadixSearch(true);
                for (int bits : RADIX_BITS_VALUES) {
                    store.setRadixBits(bits);
                    for (long q : queries) store.getLong(q, posStore);
                    for (long q : queries) store.getLongInlined(q, posStore);
                }
            }

            // --- Timed rounds ---
            long[] binaryNs = new long[BENCH_ROUNDS];
            long[][] radixNs = new long[nWidths][BENCH_ROUNDS];        // radix via getLong, per width
            long[][] radixInlinedNs = new long[nWidths][BENCH_ROUNDS]; // radix via getLongInlined, per width
            for (int r = 0; r < BENCH_ROUNDS; r++) {
                store.setRadixSearch(false);
                long t0 = System.nanoTime();
                for (long q : queries) store.getLong(q, posStore);
                binaryNs[r] = System.nanoTime() - t0;

                store.setRadixSearch(true);
                for (int bi = 0; bi < nWidths; bi++) {
                    store.setRadixBits(RADIX_BITS_VALUES[bi]);

                    long t1 = System.nanoTime();
                    for (long q : queries) store.getLong(q, posStore);
                    radixNs[bi][r] = System.nanoTime() - t1;

                    long t2 = System.nanoTime();
                    for (long q : queries) store.getLongInlined(q, posStore);
                    radixInlinedNs[bi][r] = System.nanoTime() - t2;
                }
            }

            double binMs = trimmedMeanMs(binaryNs);

            System.out.printf("%n=== Radix-width sweep vs Binary Search (k=%d) ===%n", k);
            System.out.printf("  Store: k=%d, enforceLarge=true, %,d k-mers (~%.1f MB k-mer data)%n",
                    k, entries, entries * 8.0 / (1024 * 1024));
            System.out.printf("  log2(%,d) = %.1f probes (binary)%n", entries, Math.log(entries) / Math.log(2));
            System.out.printf("  Queries per round: %,d (50%% hits, 50%% misses)%n", QUERIES_PER_ROUND);
            System.out.printf("  Rounds: %d timed, min+max discarded, %d averaged%n",
                    BENCH_ROUNDS, BENCH_ROUNDS - 2);
            System.out.printf("  Binary search:               %7.3f ms  (%5.1f ns/query)%n",
                    binMs, binMs * 1_000_000 / QUERIES_PER_ROUND);
            for (int bi = 0; bi < nWidths; bi++) {
                int bits = RADIX_BITS_VALUES[bi];
                double radMs    = trimmedMeanMs(radixNs[bi]);
                double radInlMs = trimmedMeanMs(radixInlinedNs[bi]);
                double speedup        = binMs / radMs;
                double speedupInlined = binMs / radInlMs;
                System.out.printf("  Radix %2d-bit:                %7.3f ms  (%5.1f ns/query)  %.3fx%n",
                        bits, radMs, radMs * 1_000_000 / QUERIES_PER_ROUND, speedup);
                System.out.printf("  Radix %2d-bit + inlined:      %7.3f ms  (%5.1f ns/query)  %.3fx%n",
                        bits, radInlMs, radInlMs * 1_000_000 / QUERIES_PER_ROUND, speedupInlined);

                // Sanity bound only — speedup direction depends on available cache.
                assertTrue("Radix speedup out of plausible range (k=" + k + ", bits=" + bits + "): " + speedup,
                        speedup >= 0.2 && speedup <= 5.0);
                assertTrue("Radix+inlined speedup out of plausible range (k=" + k + ", bits=" + bits + "): " + speedupInlined,
                        speedupInlined >= 0.2 && speedupInlined <= 5.0);
            }
        }
    }

    /** Mean of the timings in ms after discarding the fastest and slowest round. */
    private static double trimmedMeanMs(long[] ns) {
        long[] sorted = ns.clone();
        Arrays.sort(sorted);
        long sum = 0;
        int effective = sorted.length - 2;
        for (int r = 1; r <= effective; r++) {
            sum += sorted[r];
        }
        return sum / (1_000_000.0 * effective);
    }

    // Avoid retesting of stuff from super class.
    @Test
    public void testIdenticalResults() {
    }

    // Avoid retesting of stuff from super class.
    @Test
    public void benchmarkSpeedup() {
    }
}