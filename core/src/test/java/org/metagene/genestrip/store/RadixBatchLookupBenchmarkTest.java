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
package org.metagene.genestrip.store;

import java.io.BufferedReader;
import java.io.File;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.atomic.AtomicReference;

import org.junit.Assume;
import org.junit.Test;
import org.metagene.genestrip.io.StreamProvider;
import org.metagene.genestrip.store.RadixKMerStore.BatchBuffers;
import org.metagene.genestrip.tax.SmallTaxTree.SmallTaxIdNode;
import org.metagene.genestrip.util.CGAT;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

/**
 * Measures what {@link RadixKMerStore#getBatch} is worth over the per-k-mer
 * {@link RadixKMerStore#getLong} that {@code FastqKMerMatcher} currently uses, before the matcher's
 * hot loop is restructured to feed it.
 * <p>
 * A lookup is dominated by cache-missing loads (the radix table, the pre-filter, then the bucket's
 * binary-search steps). Issued one k-mer at a time they serialize: each miss waits for the previous.
 * {@code getBatch} runs the same work in passes across a whole batch, so many independent misses are
 * in flight at once - memory-level parallelism the per-k-mer path cannot get. The reads of a FASTQ
 * file offer this for free: a 301 bp read holds ~271 k-mers whose lookups are entirely independent.
 * <p>
 * The benchmark therefore replays <em>real</em> k-mers from a FASTQ file against a <em>real</em>
 * database and compares:
 * <ol>
 * <li>{@code getLong} - one k-mer at a time, i.e. today's matcher.</li>
 * <li>{@code getBatch(n)} - the same k-mers in batches of n, for several n.</li>
 * </ol>
 * It also asserts both paths find exactly the same k-mers, so a throughput difference cannot come
 * from one path doing less work.
 * <p>
 * Requires a database; it is skipped when none is given:
 *
 * <pre>
 * -Dgenestrip.bench.db=/path/to/viral_db.zip     the database to look up in (required)
 * -Dgenestrip.bench.fastq=/path/to/reads.fastq.gz  reads to take k-mers from (required)
 * -Dgenestrip.bench.batchKMers=20000000          how many k-mers to replay per thread and variant
 * -Dgenestrip.bench.batchSizes=32,128,512        batch sizes to try
 * -Dgenestrip.bench.batchThreads=1,9             thread counts to measure
 * </pre>
 */
public class RadixBatchLookupBenchmarkTest {
    private static final String DB_PROP = "genestrip.bench.db";
    private static final String FASTQ_PROP = "genestrip.bench.fastq";
    private static final int KMERS = Integer.getInteger("genestrip.bench.batchKMers", 20_000_000);
    private static final String BATCH_SIZES = System.getProperty("genestrip.bench.batchSizes", "32,128,512");
    private static final String THREADS = System.getProperty("genestrip.bench.batchThreads", "1,9");

    private KMerStore<SmallTaxIdNode> store;
    private RadixKMerStore<SmallTaxIdNode> radixStore;
    private long[] kmers;

    @Test
    public void testBatchedVersusSingleLookup() throws Exception {
        String dbPath = System.getProperty(DB_PROP);
        String fastqPath = System.getProperty(FASTQ_PROP);
        Assume.assumeTrue("Pass -D" + DB_PROP + " and -D" + FASTQ_PROP + " to run this benchmark.",
                dbPath != null && fastqPath != null && new File(dbPath).exists() && new File(fastqPath).exists());

        load(new File(dbPath), new File(fastqPath));

        // Both paths must find the same k-mers - otherwise a throughput difference is meaningless.
        long singleHits = countHitsSingle(Math.min(kmers.length, 200_000));
        long batchHits = countHitsBatched(Math.min(kmers.length, 200_000), 128);
        assertEquals("batched and per-k-mer lookup must find the same k-mers", singleHits, batchHits);
        assertTrue("the sample should contain hits, otherwise the benchmark measures only misses",
                singleHits > 0);
        System.out.println("  hit rate: " + (100 * singleHits / Math.min(kmers.length, 200_000)) + "%");
        System.out.println();

        int[] batchSizes = parseInts(BATCH_SIZES);
        System.out.println(String.format("%-16s %8s %12s %14s %12s", "variant", "threads", "seconds",
                "lookups/s", "vs getLong"));
        for (int threads : parseInts(THREADS)) {
            double singleSeconds = runSingle(threads, 1);
            singleSeconds = runSingle(threads, KMERS);
            System.out.println(String.format("%-16s %8d %12.2f %14.0f %11.2fx", "getLong", threads, singleSeconds,
                    (double) KMERS * threads / singleSeconds, 1.0));
            for (int batchSize : batchSizes) {
                runBatched(threads, batchSize, batchSize * 4);
                double seconds = runBatched(threads, batchSize, KMERS);
                System.out.println(String.format("%-16s %8d %12.2f %14.0f %11.2fx", "getBatch(" + batchSize + ")",
                        threads, seconds, (double) KMERS * threads / seconds, singleSeconds / seconds));
            }
            System.out.println();
        }
    }

    // --- The two lookup paths ------------------------------------------------

    private double runSingle(int threads, int kmersPerThread) throws Exception {
        return run(threads, (offset, count, sink) -> {
            long[] posStore = new long[1];
            long hits = 0;
            for (int i = 0; i < count; i++) {
                if (store.getLong(kmers[(offset + i) % kmers.length], posStore) != null) {
                    hits += posStore[0];
                }
            }
            sink[0] = hits;
        }, kmersPerThread);
    }

    private double runBatched(int threads, int batchSize, int kmersPerThread) throws Exception {
        return run(threads, (offset, count, sink) -> {
            BatchBuffers buffers = new BatchBuffers(batchSize);
            long[] hits = new long[1];
            for (int i = 0; i < count; i++) {
                if (buffers.add(kmers[(offset + i) % kmers.length], i)) {
                    radixStore.getBatch(buffers, (kmer, payload, value) -> hits[0] += payload);
                }
            }
            // The trailing partial batch still has to be looked up.
            if (!buffers.isEmpty()) {
                radixStore.getBatch(buffers, (kmer, payload, value) -> hits[0] += payload);
            }
            sink[0] = hits[0];
        }, kmersPerThread);
    }

    // --- Harness -------------------------------------------------------------

    private interface Work {
        void run(int offset, int count, long[] sink);
    }

    private double run(int threads, Work work, int kmersPerThread) throws Exception {
        CountDownLatch start = new CountDownLatch(1);
        CountDownLatch done = new CountDownLatch(threads);
        AtomicReference<Throwable> failure = new AtomicReference<>();
        long[] sinks = new long[threads * 16];
        Thread[] ts = new Thread[threads];
        for (int t = 0; t < threads; t++) {
            final int index = t;
            ts[t] = new Thread(() -> {
                try {
                    long[] sink = new long[1];
                    int offset = (int) ((long) index * kmers.length / Math.max(1, threads));
                    start.await();
                    work.run(offset, kmersPerThread, sink);
                    sinks[index * 16] = sink[0];
                } catch (Throwable th) {
                    failure.compareAndSet(null, th);
                } finally {
                    done.countDown();
                }
            });
        }
        for (Thread t : ts) {
            t.start();
        }
        long t0 = System.nanoTime();
        start.countDown();
        done.await();
        double seconds = (System.nanoTime() - t0) / 1e9;
        for (Thread t : ts) {
            t.join();
        }
        if (failure.get() != null) {
            throw new RuntimeException(failure.get());
        }
        return seconds;
    }

    // --- Fixture -------------------------------------------------------------

    @SuppressWarnings("unchecked")
    private void load(File dbFile, File fastq) throws Exception {
        System.out.println("Radix batched vs per-k-mer lookup");
        System.out.println("  database: " + dbFile);
        long t0 = System.currentTimeMillis();
        Database database = Database.load(dbFile, true);
        store = database.convertKMerStore();
        store.setUseFilter(true);
        Assume.assumeTrue("This benchmark needs a RadixKMerStore; this database holds a "
                + store.getClass().getSimpleName(), store instanceof RadixKMerStore);
        radixStore = (RadixKMerStore<SmallTaxIdNode>) store;
        System.out.println("  loaded in " + ((System.currentTimeMillis() - t0) / 1000) + " s, "
                + store.getEntries() + " entries, k=" + store.getK());

        int k = store.getK();
        List<Long> collected = new ArrayList<>();
        int[] badPos = new int[1];
        try (InputStream in = StreamProvider.getInputStreamForFile(fastq);
             BufferedReader br = new BufferedReader(new InputStreamReader(in, StandardCharsets.UTF_8))) {
            int lineNo = 0;
            for (String line = br.readLine(); line != null && collected.size() < KMERS; line = br.readLine()) {
                if (lineNo++ % 4 != 1) {
                    continue;
                }
                byte[] read = line.getBytes(StandardCharsets.US_ASCII);
                // Same k-mer walk the matcher does, so the hit/miss mix is the real one.
                for (int i = 0; i + k <= read.length && collected.size() < KMERS; i++) {
                    long straight = CGAT.kMerToLongStraight(read, i, k, badPos);
                    if (straight == -1) {
                        i = badPos[0];
                        continue;
                    }
                    long reverse = CGAT.kMerToLongReverse(read, i, k, null);
                    collected.add(CGAT.standardKMer(straight, reverse));
                }
            }
        }
        kmers = new long[collected.size()];
        for (int i = 0; i < kmers.length; i++) {
            kmers[i] = collected.get(i);
        }
        System.out.println("  k-mers:   " + kmers.length + " taken from " + fastq.getName());
    }

    private long countHitsSingle(int count) {
        long[] posStore = new long[1];
        long hits = 0;
        for (int i = 0; i < count; i++) {
            if (store.getLong(kmers[i], posStore) != null) {
                hits++;
            }
        }
        return hits;
    }

    private long countHitsBatched(int count, int batchSize) {
        BatchBuffers buffers = new BatchBuffers(batchSize);
        long[] hits = new long[1];
        for (int i = 0; i < count; i++) {
            if (buffers.add(kmers[i], 1)) {
                radixStore.getBatch(buffers, (kmer, payload, value) -> hits[0]++);
            }
        }
        if (!buffers.isEmpty()) {
            radixStore.getBatch(buffers, (kmer, payload, value) -> hits[0]++);
        }
        return hits[0];
    }

    private int[] parseInts(String value) {
        String[] parts = value.split(",");
        int[] result = new int[parts.length];
        for (int i = 0; i < parts.length; i++) {
            result[i] = Integer.parseInt(parts[i].trim());
        }
        return result;
    }
}
