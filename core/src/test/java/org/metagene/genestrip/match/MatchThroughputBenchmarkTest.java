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

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.lang.management.ManagementFactory;
import java.lang.management.OperatingSystemMXBean;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.logging.impl.SimpleLog;
import org.junit.Assume;
import org.junit.Test;
import org.metagene.genestrip.APITest;
import org.metagene.genestrip.DefaultExecutionContext;
import org.metagene.genestrip.ExecutionContext;
import org.metagene.genestrip.GSCommon;
import org.metagene.genestrip.GSConfigKey;
import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.fastq.AbstractLoggingFastqStreamer;
import org.metagene.genestrip.io.StreamProvider;
import org.metagene.genestrip.io.StreamingFileResource;
import org.metagene.genestrip.io.StreamingResourceListStream;
import org.metagene.genestrip.io.StreamingResourceStream;
import org.metagene.genestrip.store.Database;
import org.metagene.genestrip.store.KMerStore;
import org.metagene.genestrip.tax.SmallTaxTree;
import org.metagene.genestrip.tax.SmallTaxTree.SmallTaxIdNode;
import org.metagene.genestrip.util.GSLogFactory;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

/**
 * Locates the throughput ceiling of the classification pipeline by timing the very same input
 * through progressively more of it. It answers the question why the CPUs do not saturate when the
 * matcher runs with many consumer threads: everything before the consumers - reading the file,
 * inflating the GZIP stream and parsing the FASTQ records - happens on the single producer thread,
 * so the consumers can never run faster than that thread feeds them.
 * <p>
 * The stages, each streaming the identical bytes:
 * <ol>
 * <li>{@code raw read} - the compressed file is read block-wise and discarded, no inflating, no
 *     parsing. This is the disk (or page cache) ceiling.</li>
 * <li>{@code gunzip} - the same file through {@link java.util.zip.GZIPInputStream}, bytes
 *     discarded. The drop against {@code raw read} is the cost of inflating.</li>
 * <li>{@code parse} - a no-op {@link AbstractLoggingFastqStreamer} that parses every
 *     FASTQ record but does nothing with it, run inline without consumer threads. The drop against
 *     {@code gunzip} is the cost of the record parser.</li>
 * <li>{@code parse (n threads)} - the same no-op streamer, but handing every read to {@code n}
 *     consumer threads that immediately drop it. The difference to the previous stage is the pure
 *     queueing/handoff overhead.</li>
 * <li>{@code match (n threads)} - the real {@link FastqKMerMatcher} against the project database
 *     (viral by default), writing no output files. What it adds over {@code parse (n threads)} is
 *     the actual k-mer lookup and read classification work.</li>
 * </ol>
 * The report prints, per stage, the wall time, the throughput in compressed input MB/s, the reads
 * per second and - the decisive number - the average number of busy cores (process CPU time divided
 * by wall time). If {@code match} runs at nearly the same MB/s as {@code parse} while using only a
 * couple of busy cores, the producer thread, not the consumers and not the disk, is the limit.
 * <p>
 * Run it with (all properties optional):
 *
 * <pre>
 * -Dgenestrip.bench.fastq=/path/to/big.fastq.gz  input; defaults to the bundled human_virus sample
 * -Dgenestrip.bench.project=viral                project supplying the database and the match config
 * -Dgenestrip.bench.baseDir=/path/to/data        base dir holding that project; default: bundled samples
 * -Dgenestrip.bench.db=/path/to/db.zip           database file; defaults to the project's db file
 * -Dgenestrip.bench.threads=-1                   consumer threads, -1 = availableProcessors() - 1
 * -Dgenestrip.bench.threadSweep=1,2,4,8,19       repeat the match stage at each of these thread counts
 * -Dgenestrip.bench.targetBytes=250000000        compressed bytes to stream (input repeated as needed)
 * -Dgenestrip.bench.repeats=n                    repeat the input exactly n times (overrides targetBytes)
 * -Dgenestrip.bench.dir=/path/to/dir             where the repeated input is written; default: build dir
 * -Dgenestrip.bench.stages=raw,gunzip,parse,parsemt,match   subset of stages to run
 * </pre>
 *
 * A thread sweep answers directly where matching stops scaling; it loads the database once and runs
 * every point against it.
 * <p>
 * If the input is smaller than the target size, it is repeated - concatenated into one file, since
 * GZIP members concatenate - so that the measurement is not dominated by JIT warm-up and thread
 * startup. Give it a real, large FASTQ file to get realistic matching numbers: repeated data hits
 * the same database entries over and over and therefore flatters the lookup caches.
 * <p>
 * Note on the disk: unless the input is clearly larger than the machine's free RAM, the first stage
 * pulls it into the page cache and the later stages read from memory. To measure actual disk
 * behaviour, point {@code genestrip.bench.fastq} at a file bigger than RAM, or run a single stage at
 * a time via {@code genestrip.bench.stages} on a cold cache.
 */
public class MatchThroughputBenchmarkTest {
    private static final String PROP_PREFIX = "genestrip.bench.";

    // Input fastq (.gz); defaults to the small sample shipped with the release, streamed repeatedly.
    private static final String FASTQ_PROP = PROP_PREFIX + "fastq";
    // Project supplying both the database and the matcher's configuration.
    private static final String PROJECT_PROP = PROP_PREFIX + "project";
    // Base directory holding the projects; defaults to the one of the bundled sample projects.
    private static final String BASE_DIR_PROP = PROP_PREFIX + "baseDir";
    // Database file; defaults to the project's db file.
    private static final String DB_PROP = PROP_PREFIX + "db";
    // Stages to run, comma separated. Default: all of them.
    private static final String STAGES_PROP = PROP_PREFIX + "stages";
    // Number of consumer threads for the multi-threaded stages. -1 = availableProcessors() - 1.
    private static final String THREADS_PROP = PROP_PREFIX + "threads";
    // How many compressed bytes to stream per stage (the input is repeated to reach this).
    private static final String TARGET_BYTES_PROP = PROP_PREFIX + "targetBytes";
    // Explicit number of repetitions of the input, overriding TARGET_BYTES_PROP.
    private static final String REPEATS_PROP = PROP_PREFIX + "repeats";
    // Directory for the generated (repeated) input file; defaults to the build directory.
    private static final String DIR_PROP = PROP_PREFIX + "dir";
    // Comma separated thread counts at which to repeat the match stage, e.g. "1,2,4,8,19".
    private static final String THREAD_SWEEP_PROP = PROP_PREFIX + "threadSweep";

    private static final String DEFAULT_PROJECT = "viral";
    private static final long DEFAULT_TARGET_BYTES = 250L * 1000 * 1000;

    private static final int READ_BUFFER_SIZE = 1024 * 1024;

    // Cache line padding for the per-consumer sinks of the no-op streamer, so that the consumers do
    // not fight over one cache line and thereby distort the measurement they are part of.
    private static final int SINK_STRIDE = 16;

    /** The pipeline stages that can be timed, in the order in which they run. */
    private enum Stage {
        RAW("raw"), GUNZIP("gunzip"), PARSE("parse"), PARSE_MT("parsemt"), MATCH("match");

        private final String name;

        Stage(String name) {
            this.name = name;
        }
    }

    @Test
    public void testClassificationThroughput() throws Exception {
        File original = getFastqFile();
        assertTrue("Benchmark input does not exist: " + original, original.exists());

        GSProject project = new GSProject(new GSCommon(getBaseDir()),
                System.getProperty(PROJECT_PROP, DEFAULT_PROJECT), true);
        File dbFile = getDBFile(project);

        int repeats = getRepeats(original);
        File fastq = repeats > 1 ? concatenated(original, repeats) : original;
        long compressedBytes = fastq.length();
        int configuredThreads = Integer.getInteger(THREADS_PROP, project.intConfigValue(GSConfigKey.THREADS));
        // A negative thread count means "as many as the machine allows"; it is resolved here so that
        // every stage and every printed number refers to the same, actual number of threads.
        int threads = configuredThreads < 0 ? Runtime.getRuntime().availableProcessors() - 1 : configuredThreads;
        List<Integer> matchThreads = getMatchThreads(threads);

        System.out.println("Genestrip classification throughput benchmark");
        System.out.println("  input:      " + fastq + " (" + toMB(compressedBytes) + " MB compressed)");
        if (repeats > 1) {
            System.out.println("              " + repeats + " concatenated copies of " + original);
        }
        System.out.println("  project:    " + project.getName());
        System.out.println("  database:   " + (dbFile != null && dbFile.exists() ? dbFile.toString() : "<missing>"));
        System.out.println("  threads:    " + threads + (configuredThreads < 0 ? " (from " + configuredThreads + ")" : "")
                + ", matching at " + matchThreads);
        System.out.println("  processors: " + Runtime.getRuntime().availableProcessors());
        // The matching parameters come from the project's configuration, so they are printed to keep
        // a report self-describing (and to make a deliberately changed setting visible).
        System.out.println("  config:     k=" + project.intConfigValue(GSConfigKey.KMER_SIZE)
                + ", useBloomFilterForMatch=" + project.booleanConfigValue(GSConfigKey.USE_BLOOM_FILTER_FOR_MATCH)
                + ", threadQueueSize=" + project.intConfigValue(GSConfigKey.THREAD_QUEUE_SIZE)
                + ", maxClassificationPaths=" + project.intConfigValue(GSConfigKey.MAX_CLASSIFICATION_PATHS));
        System.out.println();

        List<Result> results = new ArrayList<>();
        // The per-file logging of the streamers would produce one block of lines per repetition,
        // which is neither readable nor free, so it is silenced for the duration of the benchmark.
        String logLevel = logLevelName(GSLogFactory.getInstance().getLogLevel());
        GSLogFactory.getInstance().setLogLevel("warn");
        try {
            if (isStageEnabled(Stage.RAW)) {
                results.add(runRawRead(fastq));
            }
            if (isStageEnabled(Stage.GUNZIP)) {
                results.add(runGunzip(fastq));
            }
            if (isStageEnabled(Stage.PARSE)) {
                results.add(runNoop(fastq, project, 0));
            }
            if (isStageEnabled(Stage.PARSE_MT)) {
                results.add(runNoop(fastq, project, threads));
            }
            if (isStageEnabled(Stage.MATCH)) {
                // An assumption failure alone would leave the run looking green and silent, so the
                // reason is stated explicitly before the test bows out.
                if (dbFile == null || !dbFile.exists()) {
                    System.out.println("No database at " + dbFile + " - build it or pass -D" + DB_PROP
                            + "=<db.zip>. Skipping the match stage; the read stages are reported below.");
                }
                Assume.assumeTrue("No database at " + dbFile, dbFile != null && dbFile.exists());
                // The tax tree is read-only while matching (each matcher keeps its own per-consumer
                // vote counters), so one loaded database serves every point of the sweep.
                Database database = loadDatabase(dbFile, project);
                for (int n : matchThreads) {
                    results.add(runMatch(fastq, project, database, n));
                }
            }
        } finally {
            GSLogFactory.getInstance().setLogLevel(logLevel);
            printReport(results);
        }

        // Every stage must have consumed the whole input, and every stage that parses reads must
        // have seen the same reads - otherwise the stages are not comparable and the report lies.
        long reads = -1;
        for (Result result : results) {
            assertEquals("Stage '" + result.name + "' did not consume the whole input.", compressedBytes,
                    result.inBytes);
            if (result.reads >= 0) {
                if (reads >= 0) {
                    assertEquals("Stage '" + result.name + "' saw a different number of reads.", reads, result.reads);
                }
                reads = result.reads;
            }
        }
    }

    // --- Stages --------------------------------------------------------------

    /**
     * Reads the compressed bytes and discards them: no inflating, no parsing. This is the ceiling
     * imposed by the storage layer (or by the page cache, if the input fits into it).
     */
    private Result runRawRead(File fastq) throws IOException {
        byte[] buffer = new byte[READ_BUFFER_SIZE];
        long inBytes = 0;
        long sink = 0;

        Watch watch = Watch.start();
        try (InputStream in = new FileInputStream(fastq)) {
            for (int c = in.read(buffer); c > 0; c = in.read(buffer)) {
                inBytes += c;
                sink += buffer[c - 1];
            }
        }
        return watch.stop("raw read", 0, inBytes, inBytes, -1, sink);
    }

    /**
     * Inflates the input and discards the decompressed bytes. Against {@code raw read} this is the
     * price of GZIP decompression, which runs on the producer thread only.
     */
    private Result runGunzip(File fastq) throws IOException {
        byte[] buffer = new byte[READ_BUFFER_SIZE];
        long outBytes = 0;
        long sink = 0;

        Watch watch = Watch.start();
        try (InputStream in = StreamProvider.getInputStreamForFile(fastq)) {
            for (int c = in.read(buffer); c > 0; c = in.read(buffer)) {
                outBytes += c;
                sink += buffer[c - 1];
            }
        }
        return watch.stop("gunzip", 0, fastq.length(), outBytes, -1, sink);
    }

    /**
     * Runs the full reader stack - inflating plus FASTQ record parsing - with a consumer that does
     * nothing, so the measured time is what the producer thread costs before any matching happens.
     *
     * @param threads the number of consumer threads; {@code 0} processes every read inline on the
     *                producer thread, so the difference to a positive value is the queueing overhead
     */
    private Result runNoop(File fastq, GSProject project, int threads) throws IOException {
        ExecutionContext bundle = newExecutionContext(project, threads);
        NoopFastqStreamer streamer = new NoopFastqStreamer(kOf(project), project, bundle);
        try {
            Watch watch = Watch.start();
            streamer.processFastqStreams(streamOf(fastq));
            return watch.stop("parse", threads, fastq.length(), streamer.getTotalBPs(), streamer.getTotalReads(),
                    streamer.getSink());
        } finally {
            streamer.dump();
            bundle.dump();
        }
    }

    /** Loads the database and reports how long that took - it is not part of any measured stage. */
    private Database loadDatabase(File dbFile, GSProject project) throws IOException, ClassNotFoundException {
        System.out.println("Loading database " + dbFile + " ...");
        long start = System.currentTimeMillis();
        Database database = Database.load(dbFile, project.booleanConfigValue(GSConfigKey.USE_BLOOM_FILTER_FOR_MATCH));
        System.out.println("Loaded in " + ((System.currentTimeMillis() - start) / 1000) + " s.");
        return database;
    }

    /**
     * Runs the real matcher over the same input, configured exactly as the {@code match} goal would
     * configure it, but without writing any output file, so that only the classification work itself
     * is added on top of the previous stage.
     */
    private Result runMatch(File fastq, GSProject project, Database database, int threads) throws IOException {
        boolean useFilter = project.booleanConfigValue(GSConfigKey.USE_BLOOM_FILTER_FOR_MATCH);
        SmallTaxTree taxTree = database.getTaxTree();
        KMerStore<SmallTaxIdNode> store = database.convertKMerStore();
        store.setUseFilter(useFilter);

        ExecutionContext bundle = newExecutionContext(project, threads);
        FastqKMerMatcher matcher = newMatcher(store, taxTree, project, bundle, database);
        try {
            Watch watch = Watch.start();
            MatchingResult result = matcher.runMatcher(streamOf(fastq), null, null);
            CountsPerTaxid stats = result.getGlobalStats();
            return watch.stop("match", threads, fastq.length(), stats.getReadsBPs(), stats.getReads(),
                    stats.getKMers());
        } finally {
            matcher.dump();
            bundle.dump();
        }
    }

    // --- Wiring --------------------------------------------------------------

    private StreamingResourceStream streamOf(File fastq) {
        return new StreamingResourceListStream(new StreamingFileResource(fastq));
    }

    private ExecutionContext newExecutionContext(GSProject project, int threads) {
        return new DefaultExecutionContext(null, threads, project.longConfigValue(GSConfigKey.LOG_PROGRESS_UPDATE_CYCLE));
    }

    private FastqKMerMatcher newMatcher(KMerStore<SmallTaxIdNode> store, SmallTaxTree taxTree, GSProject project,
                                        ExecutionContext bundle, Database database) {
        return new FastqKMerMatcher(store, project.intConfigValue(GSConfigKey.INITIAL_READ_SIZE_BYTES),
                project.intConfigValue(GSConfigKey.THREAD_QUEUE_SIZE), bundle,
                project.booleanConfigValue(GSConfigKey.WITH_PROBS), taxTree,
                project.intConfigValue(GSConfigKey.MAX_CLASSIFICATION_PATHS),
                project.doubleConfigValue(GSConfigKey.MAX_READ_TAX_ERROR_COUNT),
                project.doubleConfigValue(GSConfigKey.MAX_READ_CLASS_ERROR_COUNT),
                project.booleanConfigValue(GSConfigKey.WRITE_ALL),
                project.intConfigValue(GSConfigKey.MIN_KMERS_FOR_CLASS),
                database.getConfigInfo().getProperty(GSProject.DB_MD5)) {
            @Override
            protected boolean isOwnUniqueKMerBits() {
                return project.booleanConfigValue(GSConfigKey.PARALLEL_DB_MATCHING);
            }

            @Override
            protected boolean isProgressBar() {
                return false;
            }
        };
    }

    /**
     * An {@link AbstractLoggingFastqStreamer} that parses every read but does nothing with it. Only
     * the last byte of each read is accumulated into a per-consumer sink, so that neither the JIT
     * compiler can drop the parsing nor the consumers contend on a shared counter.
     */
    private static class NoopFastqStreamer extends AbstractLoggingFastqStreamer {
        private final long[] sinks;

        public NoopFastqStreamer(int k, GSProject project, ExecutionContext bundle) {
            super(k, project.intConfigValue(GSConfigKey.INITIAL_READ_SIZE_BYTES),
                    project.intConfigValue(GSConfigKey.THREAD_QUEUE_SIZE), bundle,
                    project.booleanConfigValue(GSConfigKey.WITH_PROBS));
            sinks = new long[Math.max(1, bundle.getThreads()) * SINK_STRIDE];
        }

        @Override
        protected void nextEntry(ReadEntry readStruct, ConsumerRunnable consumer) {
            if (readStruct.readSize > 0) {
                sinks[consumer.getIndex() * SINK_STRIDE] += readStruct.read[readStruct.readSize - 1];
            }
        }

        @Override
        protected boolean isProgressBar() {
            return false;
        }

        public long getTotalReads() {
            return totalReads;
        }

        public long getTotalBPs() {
            return totalBPs;
        }

        public long getSink() {
            long sum = 0;
            for (long sink : sinks) {
                sum += sink;
            }
            return sum;
        }
    }

    // --- Configuration -------------------------------------------------------

    // The base directory holding the projects; point it at another installation (e.g. the experiment
    // harness) to benchmark with that installation's project configuration and databases.
    private File getBaseDir() {
        String path = System.getProperty(BASE_DIR_PROP);
        return path != null ? new File(path) : APITest.getBaseDir();
    }

    private File getFastqFile() {
        String path = System.getProperty(FASTQ_PROP);
        return path != null ? new File(path)
                : new File(APITest.getBaseDir(), "projects/human_virus/fastq/sample.fastq.gz");
    }

    private File getDBFile(GSProject project) {
        String path = System.getProperty(DB_PROP);
        return path != null ? new File(path) : project.getDBFile();
    }

    /**
     * Writes (once, then reuses) a single input file holding {@code repeats} copies of the given
     * file. GZIP members simply concatenate, so this is a plain byte-wise copy and the result reads
     * back as one continuous FASTQ stream.
     * <p>
     * The copies are concatenated into <em>one</em> file rather than streamed as a list of resources
     * on purpose: {@link org.metagene.genestrip.fastq.AbstractFastqReader#readFastq} drains its read
     * pool by polling with {@code Thread.sleep(100)} once per file, so a list of n small files would
     * add up to n * 100 ms of idle time to every multi-threaded stage and hide what is measured
     * here. (For real inputs of a handful of large files that per-file cost is irrelevant.)
     */
    private File concatenated(File fastq, int repeats) throws IOException {
        File dir = getBenchDir();
        File target = new File(dir, fastq.getName().replaceFirst("(\\.fastq)", "-x" + repeats + "$1"));
        long expectedSize = fastq.length() * (long) repeats;
        if (target.exists() && target.length() == expectedSize) {
            return target;
        }
        System.out.println("Writing benchmark input " + target + " (" + toMB(expectedSize) + " MB) ...");
        byte[] buffer = new byte[READ_BUFFER_SIZE];
        try (OutputStream out = new FileOutputStream(target)) {
            for (int i = 0; i < repeats; i++) {
                try (InputStream in = new FileInputStream(fastq)) {
                    for (int c = in.read(buffer); c > 0; c = in.read(buffer)) {
                        out.write(buffer, 0, c);
                    }
                }
            }
        }
        return target;
    }

    private File getBenchDir() throws IOException {
        String buildDir = System.getProperty("buildDirectory", System.getProperty("java.io.tmpdir"));
        File dir = new File(System.getProperty(DIR_PROP, new File(buildDir, "bench").toString()));
        if (!dir.exists() && !dir.mkdirs()) {
            throw new IOException("Could not create benchmark directory " + dir);
        }
        return dir;
    }

    // The input is repeated as often as needed to cover the target number of compressed bytes, so
    // that the small sample file shipped with the release still yields a measurement that is not
    // dominated by JIT warm-up and thread startup.
    private int getRepeats(File fastq) {
        Integer repeats = Integer.getInteger(REPEATS_PROP);
        if (repeats != null) {
            return Math.max(1, repeats);
        }
        long targetBytes = Long.getLong(TARGET_BYTES_PROP, DEFAULT_TARGET_BYTES);
        long size = Math.max(1, fastq.length());
        return (int) Math.max(1, Math.min(Integer.MAX_VALUE, (targetBytes + size - 1) / size));
    }

    // The thread counts at which the match stage is repeated. Sweeping them shows where the
    // classification stops scaling - which is the whole point when the CPUs no longer max out.
    private List<Integer> getMatchThreads(int threads) {
        List<Integer> result = new ArrayList<>();
        String sweep = System.getProperty(THREAD_SWEEP_PROP);
        if (sweep == null) {
            result.add(threads);
            return result;
        }
        for (String value : sweep.split(",")) {
            result.add(Math.max(0, Integer.parseInt(value.trim())));
        }
        return result;
    }

    private boolean isStageEnabled(Stage stage) {
        String stages = System.getProperty(STAGES_PROP);
        if (stages == null) {
            return true;
        }
        for (String name : stages.split(",")) {
            if (stage.name.equalsIgnoreCase(name.trim())) {
                return true;
            }
        }
        return false;
    }

    // The k of the database only affects the k-mer counter of the no-op streamer, so the configured
    // value is good enough here and saves loading the database for the read-only stages.
    private int kOf(GSProject project) {
        return project.intConfigValue(GSConfigKey.KMER_SIZE);
    }

    // --- Measuring and reporting ---------------------------------------------

    /** Wall clock and process CPU time of one stage. */
    private static class Watch {
        private final long startMS;
        private final long startCpuNanos;

        private Watch() {
            // A quiescent starting point makes the CPU time attributable to the stage rather than to
            // whatever the previous stage left running.
            System.gc();
            startCpuNanos = processCpuNanos();
            startMS = System.currentTimeMillis();
        }

        static Watch start() {
            return new Watch();
        }

        Result stop(String name, int threads, long inBytes, long outBytes, long reads, long sink) {
            long wallMS = System.currentTimeMillis() - startMS;
            long cpuNanos = processCpuNanos();
            return new Result(name, threads, wallMS, cpuNanos < 0 || startCpuNanos < 0 ? -1 : cpuNanos - startCpuNanos,
                    inBytes, outBytes, reads, sink);
        }
    }

    /** The measurement of one stage. */
    private static class Result {
        /** The stage name including the thread count, as printed in the report. */
        final String name;
        /** The stage name without the thread count, used to look a stage up. */
        final String baseName;
        final int threads;
        final long wallMS;
        /** Process CPU time consumed during the stage in nanoseconds, or {@code -1} if unavailable. */
        final long cpuNanos;
        /** Compressed input bytes consumed. */
        final long inBytes;
        /** Uncompressed bytes produced (base pairs for the parsing stages). */
        final long outBytes;
        /** Reads parsed, or {@code -1} for the stages that do not parse. */
        final long reads;
        /** Kept only so that nothing the stage computed can be optimized away. */
        final long sink;

        Result(String name, int threads, long wallMS, long cpuNanos, long inBytes, long outBytes, long reads,
               long sink) {
            this.baseName = name;
            this.name = threads > 0 ? name + " (" + threads + (threads == 1 ? " thread)" : " threads)") : name;
            this.threads = threads;
            this.wallMS = wallMS;
            this.cpuNanos = cpuNanos;
            this.inBytes = inBytes;
            this.outBytes = outBytes;
            this.reads = reads;
            this.sink = sink;
        }

        double seconds() {
            return wallMS / 1000d;
        }

        /** Compressed input MB per second - the one figure comparable across all stages. */
        double inMBPerSecond() {
            return wallMS == 0 ? Double.NaN : inBytes / 1000d / wallMS;
        }

        double outMBPerSecond() {
            return wallMS == 0 ? Double.NaN : outBytes / 1000d / wallMS;
        }

        double readsPerSecond() {
            return reads < 0 || wallMS == 0 ? Double.NaN : reads * 1000d / wallMS;
        }

        /** Average number of cores kept busy: process CPU time over wall time. */
        double busyCores() {
            return cpuNanos < 0 || wallMS == 0 ? Double.NaN : cpuNanos / 1e6 / wallMS;
        }
    }

    private void printReport(List<Result> results) {
        if (results.isEmpty()) {
            return;
        }
        System.out.println();
        System.out.println(String.format("%-22s %10s %12s %12s %12s %12s", "stage", "seconds", "in MB/s", "out MB/s",
                "reads/s", "busy cores"));
        System.out.println(String.format("%-22s %10s %12s %12s %12s %12s", "----------------------", "----------",
                "------------", "------------", "------------", "------------"));
        for (Result result : results) {
            System.out.println(String.format("%-22s %10.2f %12.1f %12.1f %12.0f %12.2f", result.name, result.seconds(),
                    result.inMBPerSecond(), result.outMBPerSecond(), result.readsPerSecond(), result.busyCores()));
        }
        System.out.println();
        printInterpretation(results);
        // Referenced so that no stage's work can be considered dead code.
        long sink = 0;
        for (Result result : results) {
            sink += result.sink;
        }
        System.out.println("(checksum " + sink + ")");
    }

    // Spells out what the numbers mean for the question the benchmark exists to answer: is the
    // producer thread, the consumers or the storage the limiting factor?
    private void printInterpretation(List<Result> results) {
        Result raw = findResult(results, "raw read", null);
        Result gunzip = findResult(results, "gunzip", null);
        Result parse = findResult(results, "parse", Boolean.FALSE);
        Result parseMT = findResult(results, "parse", Boolean.TRUE);
        // With a thread sweep there are several match results; the last one uses the most threads
        // and is the one the scaling question is about.
        Result match = findLastResult(results, "match");
        int effectiveThreads = match != null ? match.threads : 0;

        if (raw != null && gunzip != null) {
            System.out.println(String.format(
                    "Inflating costs %.1fx the time of just reading the bytes (%.1f vs %.1f MB/s in).",
                    raw.inMBPerSecond() / gunzip.inMBPerSecond(), gunzip.inMBPerSecond(), raw.inMBPerSecond()));
        }
        if (gunzip != null && parse != null) {
            System.out.println(String.format(
                    "Parsing the records on top of that costs another %.0f%% (%.1f MB/s in, %.2f busy cores).",
                    100 * (gunzip.inMBPerSecond() / parse.inMBPerSecond() - 1), parse.inMBPerSecond(),
                    parse.busyCores()));
        }
        // With a sweep, the point at which matching is fastest says where more threads stop paying.
        Result best = null;
        int matchRuns = 0;
        for (Result result : results) {
            if (result.baseName.equals("match")) {
                matchRuns++;
                if (best == null || result.wallMS < best.wallMS) {
                    best = result;
                }
            }
        }
        if (matchRuns > 1) {
            System.out.println(String.format(
                    "Of the swept thread counts, matching is fastest at %d threads (%.1f s, %.0f reads/s, %.2f busy "
                            + "cores); more consumer threads than that do not pay off.",
                    best.threads, best.seconds(), best.readsPerSecond(), best.busyCores()));
        }

        Result producer = parseMT != null ? parseMT : parse;
        if (producer != null && match != null && match.wallMS > 0) {
            // What the producer thread costs is unavoidable: it is in the matching time as well. So
            // its share of the matching time says how much of the run no number of consumers can
            // ever remove.
            double producerShare = 100d * producer.wallMS / match.wallMS;
            System.out.println(String.format(
                    "Matching with %d consumer threads took %.1f s; getting the same reads to consumers that do "
                            + "nothing at all took %.1f s, i.e. %.0f%% of the matching time is the producer thread alone.",
                    effectiveThreads, match.seconds(), producer.seconds(), producerShare));
            System.out.println(String.format(
                    "Whatever the consumers do, this input cannot be fed faster than %.1f MB/s in / %.0f reads/s, and "
                            + "matching currently runs at %.1f MB/s in / %.0f reads/s.",
                    producer.inMBPerSecond(), producer.readsPerSecond(), match.inMBPerSecond(),
                    match.readsPerSecond()));
            System.out.println(String.format("Matching kept %.2f of %d cores busy on average (producer plus %d "
                    + "consumers).", match.busyCores(), effectiveThreads + 1, effectiveThreads));
            if (producerShare > 80) {
                System.out.println("=> The single producer thread (read + inflate + parse) is the bottleneck. More "
                        + "consumer threads cannot help; processing several input files in parallel, or a faster "
                        + "(parallel) decompressor, can.");
            } else if (raw != null && raw.wallMS > 0.5 * match.wallMS) {
                System.out.println("=> Just reading the raw bytes already takes a large part of the matching time: "
                        + "for this input the machine really is limited by storage.");
            } else if (match.busyCores() >= 0.8 * (effectiveThreads + 1)) {
                System.out.println("=> The consumers are saturated and do the bulk of the work; the producer is not "
                        + "the limit for this input.");
            } else {
                System.out.println(String.format("=> Neither the producer (%.0f%% of the time) nor the consumers "
                        + "(%.2f of %d cores busy) are saturated: the consumers are stalling - on memory latency of "
                        + "the k-mer lookups, on the queue, or on shared state.", producerShare, match.busyCores(),
                        effectiveThreads + 1));
            }
        }
    }

    // Finds the (single) result of the given stage, optionally restricted to the inline (0 threads)
    // or the multi-threaded variant of it.
    private Result findLastResult(List<Result> results, String baseName) {
        Result found = null;
        for (Result result : results) {
            if (result.baseName.equals(baseName)) {
                found = result;
            }
        }
        return found;
    }

    private Result findResult(List<Result> results, String baseName, Boolean multiThreaded) {
        for (Result result : results) {
            if (result.baseName.equals(baseName)
                    && (multiThreaded == null || multiThreaded.equals(result.threads > 0))) {
                return result;
            }
        }
        return null;
    }

    private static String logLevelName(int level) {
        switch (level) {
        case SimpleLog.LOG_LEVEL_TRACE:
            return "trace";
        case SimpleLog.LOG_LEVEL_DEBUG:
            return "debug";
        case SimpleLog.LOG_LEVEL_WARN:
            return "warn";
        case SimpleLog.LOG_LEVEL_ERROR:
            return "error";
        case SimpleLog.LOG_LEVEL_FATAL:
            return "fatal";
        case SimpleLog.LOG_LEVEL_OFF:
            return "off";
        default:
            return "info";
        }
    }

    private static long processCpuNanos() {
        OperatingSystemMXBean os = ManagementFactory.getOperatingSystemMXBean();
        if (os instanceof com.sun.management.OperatingSystemMXBean) {
            return ((com.sun.management.OperatingSystemMXBean) os).getProcessCpuTime();
        }
        return -1;
    }

    private static String toMB(long bytes) {
        return String.format("%.1f", bytes / 1e6);
    }
}
