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
package org.metagene.genestrip.match;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.Arrays;
import java.util.Random;

import org.junit.Test;
import org.metagene.genestrip.store.KMerSortedArray;
import org.metagene.genestrip.store.RadixKMerStore;
import org.metagene.genestrip.util.CGAT;

/**
 * Benchmarks the binary search of {@link KMerSortedArray} against the radix-indexed lookup of
 * {@link RadixKMerStore}.
 * <p>
 * Both stores are populated with the <em>identical</em> set of k-mers: a {@link KMerSortedArray} is
 * built first, then a {@link RadixKMerStore} is filled from exactly its entries (with the duplicate
 * pre-filter disabled during the bulk load, see {@link RadixKMerStore#setUseFilter(boolean)}, so the
 * reserved buckets fill losslessly). The {@link KMerSortedArray} uses {@code enforceLarge=true} to
 * exercise its BigArray code path, the regime that matters for production databases.
 * <p>
 * For {@code k=31} the top 16 bits of a k-mer (bits 48..61) feed the radix table, so lookups skip
 * straight to a small per-prefix bucket; for smaller k all k-mers share radix bucket 0 and the radix
 * store degenerates to a single sorted array — included to show the radix index does not hurt there.
 */
public class RadixKMerStoreBenchmarkTest {

    // The benchmark is repeated for each of these k-mer sizes.
    private static final int[] K_VALUES = {16, 21, 31};
    // Three taxonomy ids spread across the store.
    private static final String[] TAXIDS = {"1", "2", "3"};

    // Radix width for the RadixKMerStore under test.
    private static final int RADIX_BITS = RadixKMerStore.DEFAULT_RADIX_BITS;

    // Large store for the benchmark. Override with -Dgenestrip.bench.radixKmersPerTaxid=<n>.
    // Two stores are held simultaneously, so this defaults lower than the parent's sweep.
    private static final int BENCH_KMERS_PER_TAXID =
            Integer.getInteger("genestrip.bench.radixKmersPerTaxid", 50_000_000);

    // Small store for the cross-store correctness check (all entries verified).
    private static final int CORRECTNESS_KMERS_PER_TAXID = 20_000; // 60 K total

    private static final int QUERIES_PER_ROUND = 100_000;

    // Built pair of stores with identical content.
    private static final class Stores {
        final KMerSortedArray<String> sorted;
        final RadixKMerStore<String> radix;
        final long entries;

        Stores(KMerSortedArray<String> sorted, RadixKMerStore<String> radix, long entries) {
            this.sorted = sorted;
            this.radix = radix;
            this.entries = entries;
        }
    }

    // --- Store builders ------------------------------------------------------

    /**
     * Builds a {@link KMerSortedArray} of random k-mers and a {@link RadixKMerStore} holding exactly
     * the same entries. {@code enforceLarge} only affects the {@link KMerSortedArray}; the radix
     * store always uses its per-bucket {@code long[]} layout.
     */
    private static Stores buildStores(int k, int kmersPerTaxid, boolean enforceLarge, Random rng) {
        long total = (long) kmersPerTaxid * TAXIDS.length;
        KMerSortedArray<String> sorted = new KMerSortedArray<>(
                k, 0.001, 0.001, Arrays.asList(TAXIDS), enforceLarge, true);
        sorted.initSize(total);

        byte[] seq = new byte[k];
        for (String taxid : TAXIDS) {
            int added = 0;
            while (added < kmersPerTaxid) {
                for (int j = 0; j < k; j++) {
                    seq[j] = CGAT.DECODE_TABLE[rng.nextInt(4)];
                }
                long canonical = CGAT.kMerToLong(seq, 0, k, null);
                if (canonical != -1 && sorted.putLong(canonical, taxid)) {
                    added++;
                }
            }
        }
        sorted.optimize();
        long entries = sorted.getEntries();

        // Count the distinct k-mers per radix prefix, then reserve exactly that much per bucket.
        int[] bucketSizes = new int[1 << RADIX_BITS];
        for (long pos = 0; pos < entries; pos++) {
            bucketSizes[RadixKMerStore.radixOf(sorted.getKMerAt(pos), RADIX_BITS)]++;
        }

        RadixKMerStore<String> radix = new RadixKMerStore<>(
                k, RADIX_BITS, bucketSizes, 0.001, 0.001, Arrays.asList(TAXIDS), true);
        // The source k-mers are already distinct: skip the probabilistic duplicate check so every
        // entry is stored (no Bloom false-positive drops), giving the two stores identical content.
        radix.setUseFilter(false);
        for (long pos = 0; pos < entries; pos++) {
            long kmer = sorted.getKMerAt(pos);
            String value = sorted.getValueForIndex(sorted.indexAtPosition(pos));
            radix.putLong(kmer, value);
        }
        radix.setUseFilter(true);
        radix.optimize();

        return new Stores(sorted, radix, entries);
    }

    /**
     * Samples {@code n} stored k-mers from uniformly random positions of the sorted store
     * (O(n), independent of store size). They are present in both stores.
     */
    private static long[] sampleHits(KMerSortedArray<?> sorted, int n, Random rng) {
        long entries = sorted.getEntries();
        long[] sample = new long[n];
        for (int i = 0; i < n; i++) {
            sample[i] = sorted.getKMerAt((rng.nextLong() >>> 1) % entries);
        }
        return sample;
    }

    /** Generates {@code n} random k-mers that are (with overwhelming probability) not stored. */
    private static long[] generateMisses(int k, int n, Random rng) {
        long[] misses = new long[n];
        byte[] seq = new byte[k];
        for (int i = 0; i < n; i++) {
            long canonical;
            do {
                for (int j = 0; j < k; j++) {
                    seq[j] = CGAT.DECODE_TABLE[rng.nextInt(4)];
                }
                canonical = CGAT.kMerToLong(seq, 0, k, null);
            } while (canonical == -1);
            misses[i] = canonical;
        }
        return misses;
    }

    /** Mean of the timings in ms after discarding the fastest and slowest round. */
    private static double trimmedMeanMs(long[] ns) {
        long[] s = ns.clone();
        Arrays.sort(s);
        long sum = 0;
        int effective = s.length - 2;
        for (int r = 1; r <= effective; r++) {
            sum += s[r];
        }
        return sum / (1_000_000.0 * effective);
    }

    // --- Correctness ---------------------------------------------------------

    /**
     * Verifies that {@link KMerSortedArray} (binary search) and {@link RadixKMerStore} return the
     * same value for every stored k-mer and agree on misses. The Bloom filter is disabled on both
     * so that misses actually reach the respective search code instead of being short-circuited.
     */
    @Test
    public void testRadixStoreIdenticalToBinary() {
        for (int k : K_VALUES) {
            Stores s = buildStores(k, CORRECTNESS_KMERS_PER_TAXID, true, new Random(54321L));
            s.sorted.setUseFilter(false);
            s.radix.setUseFilter(false);

            long[] posStore = new long[1];

            // Hits: every stored k-mer must resolve to the same value in both stores.
            for (long pos = 0; pos < s.entries; pos++) {
                long kmer = s.sorted.getKMerAt(pos);
                String binary = s.sorted.getLong(kmer, posStore);
                String radix = s.radix.getLong(kmer, posStore);
                assertEquals("k=" + k + " value mismatch for stored kmer at pos " + pos, binary, radix);
            }

            // Random queries (mostly misses): the two stores must agree.
            int misses = 0;
            byte[] seq = new byte[k];
            Random rng = new Random(98765L);
            for (int i = 0; i < 100_000; i++) {
                for (int j = 0; j < k; j++) {
                    seq[j] = CGAT.DECODE_TABLE[rng.nextInt(4)];
                }
                long kmer = CGAT.kMerToLong(seq, 0, k, null);
                if (kmer == -1) {
                    continue;
                }
                String binary = s.sorted.getLong(kmer, posStore);
                String radix = s.radix.getLong(kmer, posStore);
                assertEquals("k=" + k + " disagreement on random kmer " + i, binary, radix);
                if (binary == null) {
                    misses++;
                }
            }
            System.out.printf("%nCorrectness (k=%d): %,d stored k-mers + %,d random misses verified identical%n",
                    k, s.entries, misses);
        }
    }

    // --- Benchmark -----------------------------------------------------------

    /**
     * Times {@link KMerSortedArray#getLong} (binary search) against {@link RadixKMerStore#getLong}
     * on the same 50/50 hit/miss query stream, for each k in {@link #K_VALUES}.
     */
    @Test
    public void benchmarkRadixStoreVsBinary() {
        final int WARMUP_ITERS = 5;
        final int BENCH_ROUNDS = 7;
        final int HALF = QUERIES_PER_ROUND / 2;

        for (int k : K_VALUES) {
            Random rng = new Random(77777L);
            Stores s = buildStores(k, BENCH_KMERS_PER_TAXID, true, rng);

            long[] hits = sampleHits(s.sorted, HALF, rng);
            long[] misses = generateMisses(k, HALF, rng);
            long[] queries = new long[QUERIES_PER_ROUND];
            for (int i = 0; i < HALF; i++) {
                queries[2 * i] = hits[i];
                queries[2 * i + 1] = misses[i];
            }

            long[] posStore = new long[1];

            // --- Warm-up ---
            for (int w = 0; w < WARMUP_ITERS; w++) {
                for (long q : queries) s.sorted.getLong(q, posStore);
                for (long q : queries) s.radix.getLong(q, posStore);
            }

            // --- Timed rounds ---
            long[] binaryNs = new long[BENCH_ROUNDS];
            long[] radixNs = new long[BENCH_ROUNDS];
            for (int r = 0; r < BENCH_ROUNDS; r++) {
                long t0 = System.nanoTime();
                for (long q : queries) s.sorted.getLong(q, posStore);
                binaryNs[r] = System.nanoTime() - t0;

                long t1 = System.nanoTime();
                for (long q : queries) s.radix.getLong(q, posStore);
                radixNs[r] = System.nanoTime() - t1;
            }

            double binMs = trimmedMeanMs(binaryNs);
            double radMs = trimmedMeanMs(radixNs);
            double speedup = binMs / radMs;

            System.out.printf("%n=== KMerSortedArray binary search vs RadixKMerStore (k=%d) ===%n", k);
            System.out.printf("  Stores: %,d k-mers each (~%.1f MB k-mer data per store)%n",
                    s.entries, s.entries * 8.0 / (1024 * 1024));
            System.out.printf("  log2(%,d) = %.1f probes (binary, whole array)%n",
                    s.entries, Math.log(s.entries) / Math.log(2));
            System.out.printf("  Queries per round: %,d (50%% hits, 50%% misses)%n", QUERIES_PER_ROUND);
            System.out.printf("  Rounds: %d timed, min+max discarded, %d averaged%n",
                    BENCH_ROUNDS, BENCH_ROUNDS - 2);
            System.out.printf("  KMerSortedArray (binary):  %7.3f ms  (%5.1f ns/query)%n",
                    binMs, binMs * 1_000_000 / QUERIES_PER_ROUND);
            System.out.printf("  RadixKMerStore:            %7.3f ms  (%5.1f ns/query)  %.3fx%n",
                    radMs, radMs * 1_000_000 / QUERIES_PER_ROUND, speedup);

            // Sanity bound only — the direction of the speedup depends on k and available cache.
            assertTrue("RadixKMerStore speedup out of plausible range (k=" + k + "): " + speedup,
                    speedup >= 0.1 && speedup <= 10.0);
        }
    }
}
