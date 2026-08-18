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
import java.nio.charset.StandardCharsets;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

import org.metagene.genestrip.ExecutionContext;
import org.metagene.genestrip.fastq.AbstractFastqReader;
import org.metagene.genestrip.fastq.AbstractLoggingFastqStreamer;
import org.metagene.genestrip.io.StreamProvider;
import org.metagene.genestrip.io.StreamingResource;
import org.metagene.genestrip.io.StreamingResourceListStream;
import org.metagene.genestrip.io.StreamingResourceStream;
import org.metagene.genestrip.store.KMerStore;
import org.metagene.genestrip.util.LargeBitVector;
import org.metagene.genestrip.store.RadixKMerStore;
import org.metagene.genestrip.tax.SmallTaxTree;
import org.metagene.genestrip.tax.SmallTaxTree.SmallTaxIdNode;
import org.metagene.genestrip.util.ByteArrayUtil;
import org.metagene.genestrip.util.CGAT;


/**
 * Matches the k-mers of FASTQ (or FASTA) reads against the k-mer database and classifies
 * each read to a tax id in the style of Kraken/KrakenUniq. For every read it walks its
 * k-mers, votes on candidate taxonomic paths, resolves the classification via the lowest
 * common ancestor (subject to error and count thresholds) and accumulates per-tax-id
 * statistics into {@link CountsPerTaxid}. Optionally it writes the filtered reads and a
 * Kraken-style output file and counts unique k-mers. Reads are processed by several
 * consumer threads in parallel.
 */
public class FastqKMerMatcher extends AbstractLoggingFastqStreamer {
    /** Sentinel node marking an invalid (ambiguous) k-mer position within a read. */
    protected final static SmallTaxIdNode INVALID_NODE = new SmallTaxIdNode("INVALID", null, null);

    /** The k-mer store (database) mapping k-mers to their tax id nodes. */
    protected final KMerStore<SmallTaxIdNode> kmerStore;
    /** MD5 checksum identifying the database used for matching. */
    protected final String dbMD5;

    /**
     * This matcher's own record of the distinct k-mers it matched, one bit per store position, or
     * {@code null} when the marks are kept in the database instead (see {@link #isOwnUniqueKMerBits}).
     * Marking inside the database is cheaper - the lookup has just loaded that cache line - but it
     * writes to the database, so it rules out running several matchers against one loaded database.
     */
    protected LargeBitVector uniqueKmerBits;
    /** Per-store-index statistics accumulators, indexed by a tax id node's store index. */
    protected final CountsPerTaxid[] statsIndex;

    /**
     * The store as a {@link RadixKMerStore} when its batched lookup can be used, otherwise
     * {@code null}. Looking a read's k-mers up one at a time serializes their cache misses, because
     * each miss is only issued once the previous one has returned; the batched lookup resolves a whole
     * batch in passes, keeping many independent misses in flight (memory-level parallelism). A read's
     * k-mers are ideal for this: their lookups do not depend on each other.
     */
    private final RadixKMerStore<SmallTaxIdNode> batchStore;

    /** Number of k-mers looked up per batch. Beyond ~128 the gain flattens out (measured). */
    private static final int BATCH_SIZE = 128;

    /**
     * One matching consumer thread and everything that belongs to it alone: the per-read vote
     * counters, the batch buffers and the prefetched lookup results. Keeping them here rather than in
     * arrays indexed by a consumer index means a thread can only reach its own state, so no two
     * threads can write neighbouring slots of a shared array - the mistake that per-node counter
     * slots used to make.
     * <p>
     * The fields are filled in by the matcher's constructor (see {@link #initConsumers}), because a
     * consumer is created by the reader's constructor, before this class's own fields exist.
     */
    protected static class MatcherConsumer extends ConsumerRunnable {
        /** Vote counters of the current read, indexed by node position; {@code null} without classification. */
        protected int[] nodeCounts;
        /** Read key per counter: a counter only counts for the read whose number is stored here. */
        protected long[] nodeCountInitKeys;
        /** Last read number seen per store index, to count reads with at least one k-mer once. */
        protected long[] readNoPerStat;
        /** Buffers collecting a read's k-mers for the batched lookup; {@code null} without batching. */
        protected RadixKMerStore.BatchBuffers buffers;
        /** Prefetched node per k-mer start position of the current read, or {@code null} for a miss. */
        protected SmallTaxIdNode[] prefetchedNodes;
        /** Prefetched storage positions, valid where the corresponding node is not {@code null}. */
        protected long[] prefetchedPositions;

        /**
         * Creates the consumer for the given matcher and index.
         *
         * @param reader the matcher this consumer belongs to.
         * @param index the index of this consumer among the matcher's consumers.
         */
        public MatcherConsumer(AbstractFastqReader reader, int index) {
            super(reader, index);
        }
    }

    @Override
    protected ConsumerRunnable createRunnable(int rindex, Object... config) {
        return new MatcherConsumer(this, rindex);
    }

    /** Maximum number of candidate taxonomic paths tracked per read. */
    protected final int maxPaths;
    /** The taxonomy tree used for voting and lowest-common-ancestor resolution. */
    protected final SmallTaxTree taxTree;
    /** Maximum allowed tax error count per read before it is discarded. */
    protected final double maxReadTaxErrorCount;
    /** Maximum allowed classification error count per read before it is discarded. */
    protected final double maxReadClassErrorCount;
    /** Output stream for the filtered (matched) reads, or {@code null} if none. */
    protected OutputStream indexed;
    /** Minimum vote count threshold applied during classification. */
    protected final int threshold;

    // This should stay a box type for the line root.get(taxid.getTaxId(),
    // maxReadSize);
    /** Initial read size in bytes (kept boxed for reuse as a map key). */
    protected final Integer initialReadSize;

    // A PrintStream is implicitly synchronized. So we don't need to worry about
    // multi-threading when using it.
    /** Print stream for the Kraken-style per-read output, or {@code null} if none. */
    protected PrintStream out;
    /** Whether to write all reads to the Kraken-style output, not just classified ones. */
    protected final boolean writeAll;

    /** Number of parallel consumer threads processing reads. */
    protected int consumers;

    private AfterMatchCallback afterMatchCallback;

    /**
     * Creates a matcher over the given database and taxonomy.
     *
     * @param kmerStore             the k-mer store (database) to match against
     * @param initialReadSize       the initial read size in bytes
     * @param maxQueueSize          the maximum size of the read queue
     * @param bundle                the execution context providing the worker threads
     * @param withProbs             whether quality/probability information is processed
     * @param taxTree               the taxonomy tree used for voting and LCA resolution
     * @param maxPaths              the maximum number of candidate taxonomic paths per read
     * @param maxReadTaxErrorCount  the maximum allowed tax error count per read
     * @param maxReadClassErrorCount the maximum allowed classification error count per read
     * @param writeAll              whether to write all reads, not just classified ones
     * @param threshold             the minimum vote count threshold for classification
     * @param dbMD5                 the MD5 checksum identifying the database
     */
    public FastqKMerMatcher(KMerStore<SmallTaxIdNode> kmerStore, int initialReadSize, int maxQueueSize,
                            ExecutionContext bundle, boolean withProbs, SmallTaxTree taxTree, int maxPaths,
                            double maxReadTaxErrorCount, double maxReadClassErrorCount, boolean writeAll, int threshold, String dbMD5) {
        super(kmerStore.getK(), initialReadSize, maxQueueSize, bundle, withProbs, maxPaths);
        consumers = bundle.getThreads() <= 0 ? 1 : bundle.getThreads();
        this.kmerStore = kmerStore;
        this.statsIndex = new CountsPerTaxid[kmerStore.getNValues()];
        this.initialReadSize = initialReadSize;
        this.taxTree = taxTree;
        this.maxReadTaxErrorCount = maxReadTaxErrorCount;
        this.maxReadClassErrorCount = maxReadClassErrorCount;
        this.maxPaths = maxPaths;
        this.writeAll = writeAll;
        this.threshold = threshold;
        this.dbMD5 = dbMD5;
        // Only an optimized radix store can serve batched lookups; anything else keeps the per-k-mer
        // path (a store that is not optimized yet has no buckets to search).
        batchStore = kmerStore instanceof RadixKMerStore && kmerStore.isOptimized()
                ? (RadixKMerStore<SmallTaxIdNode>) kmerStore : null;
        initConsumers(kmerStore, initialReadSize);
    }

    /**
     * Whether this matcher counts its distinct matched k-mers in a bit vector of its own rather than
     * marking them in the database. Keeping them here leaves the database untouched, so several
     * matchers can share one - at the price of a bit per stored k-mer and an extra memory access per
     * matched k-mer. Off by default; {@code MatchResultGoal} overrides it from the configuration.
     *
     * @return whether to keep the unique-k-mer bits in this matcher
     */
    protected boolean isOwnUniqueKMerBits() {
        return false;
    }

    /**
     * Fills in the per-consumer state. The consumers themselves are created by the reader's
     * constructor, which runs before this class's fields exist, so their contents can only be
     * computed here. The consumer threads are already started but blocked on the queue, and the
     * queue's handoff publishes these writes to them.
     *
     * @param kmerStore the store being matched against
     * @param initialReadSize the initial read size in bytes, sizing the prefetch arrays
     */
    private void initConsumers(KMerStore<SmallTaxIdNode> kmerStore, int initialReadSize) {
        // Positions are established when the tree is built, loaded or restructured, so the counters
        // can be sized from a plain read: the tree is never written to from here, which keeps it
        // read-only for the whole of matching and lets several matchers share one database.
        int nodes = taxTree == null ? 0 : taxTree.getNodeCount();
        for (int i = 0; i < consumers; i++) {
            MatcherConsumer consumer = (MatcherConsumer) getConsumer(i);
            consumer.readNoPerStat = new long[kmerStore.getNValues()];
            if (taxTree != null) {
                consumer.nodeCounts = new int[nodes];
                consumer.nodeCountInitKeys = new long[nodes];
            }
            if (batchStore != null) {
                consumer.buffers = new RadixKMerStore.BatchBuffers(BATCH_SIZE);
                consumer.prefetchedNodes = new SmallTaxIdNode[initialReadSize];
                consumer.prefetchedPositions = new long[initialReadSize];
            }
        }
    }

    /**
     * Looks up every k-mer of the read in batches and records the results per k-mer start position,
     * so that {@link #matchRead} can consume them without stalling on one cache miss at a time. The
     * k-mer walk mirrors the one in {@link #matchRead} exactly - it is pure arithmetic on the read's
     * bytes and touches no memory that could miss, so repeating it is far cheaper than serializing
     * the lookups it feeds.
     *
     * @param entry the read whose k-mers are looked up
     * @param max   the number of k-mer start positions in the read
     * @param consumer the consumer whose private buffers and prefetch arrays are used
     */
    private void prefetchKMers(final MatcherReadEntry entry, final int max, final MatcherConsumer consumer) {
        if (consumer.prefetchedNodes.length < max) {
            // Reads may grow beyond the initial buffer size; grow with them and keep the arrays.
            consumer.prefetchedNodes = new SmallTaxIdNode[max];
            consumer.prefetchedPositions = new long[max];
        }
        final SmallTaxIdNode[] lnodes = consumer.prefetchedNodes;
        final long[] positions = consumer.prefetchedPositions;
        // A miss leaves no result behind, so stale nodes from the previous read must be cleared.
        Arrays.fill(lnodes, 0, max, null);

        final RadixKMerStore.BatchBuffers buffers = consumer.buffers;
        final LargeBitVector bits = uniqueKmerBits;
        final RadixKMerStore.BatchPositionConsumer<SmallTaxIdNode> resultSink = (kmer, payload, value, position) -> {
            lnodes[payload] = value;
            positions[payload] = position;
            if (bits != null) {
                // Set here rather than while classifying: this loop runs right after the batch was
                // resolved, so these writes overlap with each other instead of being spread out.
                bits.set(position);
            }
        };

        long kmer = -1;
        long reverseKmer = -1;
        for (int i = 0; i < max; i++) {
            if (kmer == -1) {
                kmer = CGAT.kMerToLongStraight(entry.read, i, k, entry.badPos);
                if (kmer == -1) {
                    i = entry.badPos[0];
                } else {
                    reverseKmer = CGAT.kMerToLongReverse(entry.read, i, k, null);
                }
            } else {
                final byte lastBase = entry.read[i + k - 1];
                kmer = CGAT.nextKMerStraight(kmer, lastBase, k);
                if (kmer == -1) {
                    i += k - 1;
                } else {
                    reverseKmer = CGAT.nextKMerReverse(reverseKmer, lastBase, k);
                }
            }
            if (kmer != -1 && buffers.add(CGAT.standardKMer(kmer, reverseKmer), i)) {
                batchStore.getBatch(buffers, resultSink);
            }
        }
        if (!buffers.isEmpty()) {
            batchStore.getBatch(buffers, resultSink);
        }
    }

    /**
     * Counts one vote of the current read for the given node, on this consumer's own counters.
     * {@code initKey} identifies the read, so a counter left from a previous read restarts instead of
     * accumulating.
     *
     * @param node    the node to count a vote for
     * @param consumer the consumer whose private counters are used
     * @param initKey the key identifying the current read
     */
    // Made final for potential inlining by JVM
    protected final void incCount(final SmallTaxIdNode node, final MatcherConsumer consumer, final long initKey) {
        final int pos = node.getPosition();
        final long[] initKeys = consumer.nodeCountInitKeys;
        final int[] counts = consumer.nodeCounts;
        if (initKeys[pos] == initKey) {
            counts[pos]++;
        } else {
            initKeys[pos] = initKey;
            counts[pos] = 1;
        }
    }

    /**
     * Sums this consumer's counters belonging to {@code initKey} along the path from the node to the
     * root.
     *
     * @param node    the node to start summing from
     * @param consumer the consumer whose private counters are used
     * @param initKey the key identifying the current read
     * @return the sum of the matching counts from the node to the root
     */
    // Made final for potential inlining by JVM
    protected final int sumCounts(SmallTaxIdNode node, final MatcherConsumer consumer, final long initKey) {
        final long[] initKeys = consumer.nodeCountInitKeys;
        final int[] counts = consumer.nodeCounts;
        int res = 0;
        while (node != null) {
            final int pos = node.getPosition();
            if (initKeys[pos] == initKey) {
                res += counts[pos];
            }
            node = node.getParent();
        }
        return res;
    }

    /**
     * Walks from the node to the root accumulating this consumer's counters belonging to
     * {@code initKey}, and returns the lowest node at which the running sum reaches {@code threshold}.
     *
     * @param node      the node to start summing from
     * @param consumer  the consumer whose private counters are used
     * @param initKey   the key identifying the current read
     * @param threshold the running sum to reach
     * @return the lowest node where the running sum reaches {@code threshold}, or {@code null}
     */
    // Made final for potential inlining by JVM
    protected final SmallTaxIdNode lowestNodeWhereSumAboveThreshold(SmallTaxIdNode node,
                                                                   final MatcherConsumer consumer,
                                                                   final long initKey, int threshold) {
        final long[] initKeys = consumer.nodeCountInitKeys;
        final int[] counts = consumer.nodeCounts;
        int res = 0;
        while (node != null) {
            final int pos = node.getPosition();
            if (initKeys[pos] == initKey) {
                res += counts[pos];
                if (res >= threshold) {
                    return node;
                }
            }
            node = node.getParent();
        }
        return null;
    }

    @Override
    protected ReadEntry createReadEntry(int initialReadSizeBytes, boolean withProbs, Object... config) {
        return new MatcherReadEntry(initialReadSizeBytes, withProbs, (int) config[0]);
    }

    /**
     * Convenience overload of {@link #runMatcher(StreamingResourceStream, File, File)}
     * that matches a single FASTQ resource.
     *
     * @param fastq             the FASTQ resource to process
     * @param filteredFile      optional file to which matched reads are written, or {@code null}
     * @param krakenOutStyleFile optional file for Kraken-style per-read output, or {@code null}
     * @return the aggregated matching result
     * @throws IOException if reading the FASTQ resource or writing the output files fails
     */
    public MatchingResult runMatcher(StreamingResource fastq, File filteredFile, File krakenOutStyleFile)
            throws IOException {
        return runMatcher(new StreamingResourceListStream(fastq), filteredFile, krakenOutStyleFile);
    }

    /**
     * Matches all reads of the given FASTQ streams against the database and returns the
     * aggregated result.
     *
     * @param fastqs            the FASTQ resources to process
     * @param filteredFile      optional file to which matched reads are written, or {@code null}
     * @param krakenOutStyleFile optional file for Kraken-style per-read output, or {@code null}
     * @return the aggregated matching result
     * @throws IOException if reading the FASTQ streams or writing the output files fails
     */
    public MatchingResult runMatcher(StreamingResourceStream fastqs, File filteredFile, File krakenOutStyleFile)
            throws IOException {
        try (OutputStream lindexed = filteredFile != null ? StreamProvider.getOutputStreamForFile(filteredFile) : null;
             // A PrintStream is implicitly synchronized. So we don't need to worry about
             // multi threading when using it.
             PrintStream lout = krakenOutStyleFile != null
                     ? new PrintStream(StreamProvider.getOutputStreamForFile(krakenOutStyleFile), false, StandardCharsets.UTF_8)
                     : null) {
            indexed = lindexed;
            out = lout;

            initStats();
            initUniqueCounter();
            processFastqStreams(fastqs);
        }
        out = null;
        indexed = null;

        Map<String, CountsPerTaxid> taxid2Stats = new HashMap<>();
        for (CountsPerTaxid stats : statsIndex) {
            if (stats != null) {
                taxid2Stats.put(stats.getTaxid(), stats);
            }
        }

        computeUniqueKmerCounts();

        return new MatchingResult(kmerStore.getK(), taxid2Stats, dbMD5, totalReads, totalKMers, totalBPs);
    }

    /**
     * Fills each tax id's {@code uniqueKmers} from the store's visited marks: the number of distinct
     * matched k-mers per store value index (which {@link #statsIndex} is indexed by too).
     * Package-private so tests can drive it directly.
     */
    void computeUniqueKmerCounts() {
        long[] uniquePerIndex = new long[statsIndex.length];
        if (uniqueKmerBits != null) {
            final LargeBitVector bits = uniqueKmerBits;
            kmerStore.visit((store, kmer, index, pos) -> {
                if (bits.get(pos)) {
                    uniquePerIndex[index]++;
                }
            });
        } else {
            kmerStore.countVisitedPerValueIndex(uniquePerIndex);
        }
        for (int vi = 0; vi < statsIndex.length; vi++) {
            if (statsIndex[vi] != null) {
                statsIndex[vi].uniqueKmers = uniquePerIndex[vi];
            }
        }
    }

    @Override
    protected void readFastq(InputStream inputStream, boolean fasta) throws IOException {
        // Read numbers restart at zero for every file, so counters still carrying a key from the
        // previous file would be mistaken for the current read's and must be invalidated. The
        // consumers have drained by this point, so nothing reads them concurrently.
        for (int i = 0; i < consumers; i++) {
            MatcherConsumer consumer = (MatcherConsumer) getConsumer(i);
            if (consumer.nodeCountInitKeys != null) {
                Arrays.fill(consumer.nodeCountInitKeys, -1);
            }
            Arrays.fill(consumer.readNoPerStat, -1);
        }
        super.readFastq(inputStream, fasta);
    }

    // Package private for testing purposes.
    void initStats() {
        Arrays.fill(statsIndex, null);
    }

    void initUniqueCounter() {
        if (isOwnUniqueKMerBits()) {
            // The database is left alone so that other matchers may read it at the same time. The
            // vector is indexed by store position, which is dense and ordered by radix bucket, so it
            // is touched in much the same pattern as the store itself.
            if (uniqueKmerBits == null) {
                uniqueKmerBits = new LargeBitVector(kmerStore.getEntries());
            } else {
                uniqueKmerBits.clear();
            }
            return;
        }
        uniqueKmerBits = null;
        // The lookup has just loaded the entry's cache line, so marking it there avoids the second
        // random memory access a separate bit vector needs. It writes to the database though, so it
        // can only be claimed by one matcher at a time.
        if (!kmerStore.setMarkVisited(true)) {
            if (kmerStore.isMarkVisited()) {
                throw new IllegalStateException(
                        "Another matcher is already marking this database. Set parallelDbMatching=true"
                                + " on every matching run that shares a database, so each keeps its own"
                                + " unique-k-mer bits.");
            }
            throw new IllegalStateException("The database cannot mark visited k-mers: it holds "
                    + kmerStore.getNValues() + " values, which reaches into the entry bit reserved for the"
                    + " mark. Rebuild it (optionally with a wider radixStoreBits) to match it.");
        }
        kmerStore.clearVisitedMarks();
    }

    @Override
    public void dump() {
        // Release the database's marking so a following matcher can claim it.
        if (!isOwnUniqueKMerBits()) {
            kmerStore.setMarkVisited(false);
        }
        super.dump();
    }

    @Override
    // Made final for potential inlining by JVM
    protected final void nextEntry(ReadEntry entry, ConsumerRunnable consumerRunnable) throws IOException {
        final MatcherConsumer consumer = (MatcherConsumer) consumerRunnable;
        MatcherReadEntry myEntry = (MatcherReadEntry) entry;
        myEntry.bufferPos = 0;

        myEntry.usedPaths = 0;
        myEntry.classNode = null;
        for (int i = 0; i < maxPaths; i++) {
            myEntry.readTaxIdNode[i] = null;
            myEntry.counts[i] = 0;
        }

        boolean found = matchRead(myEntry, consumer);
        afterMatch(myEntry, found);
        if (afterMatchCallback != null) {
            afterMatchCallback.afterMatch(myEntry, found);
        }
    }

    /**
     * Sets the callback invoked after each read has been matched.
     *
     * @param afterMatchCallback the callback to invoke, or {@code null} for none
     */
    public void setAfterMatchCallback(AfterMatchCallback afterMatchCallback) {
        this.afterMatchCallback = afterMatchCallback;
    }

    /**
     * Called after a read has been matched: writes the read to the filtered output if it
     * matched, and appends its Kraken-style line to the output stream when applicable.
     *
     * @param myEntry the matched read entry
     * @param found   whether the read matched at least one k-mer
     * @throws IOException if writing the read or its output line fails
     */
    protected void afterMatch(MatcherReadEntry myEntry, boolean found) throws IOException {
        if (found && indexed != null) {
            rewriteInput(myEntry, indexed);
        }
        if (out != null) {
            if (writeAll || myEntry.classNode != null) {
                synchronized (out) {
                    myEntry.writeMatchDetails(out);
                }
            }
        }
    }

    /**
     * Matches a single read against the database: walks its k-mers, tracks contigs of
     * k-mers belonging to the same tax id, updates per-tax-id statistics and unique-k-mer
     * counts, votes on candidate taxonomic paths and resolves the read's classification
     * via the lowest common ancestor subject to the error and count thresholds.
     *
     * @param entry the read together with its per-read working state
     * @param index the consumer thread index (selects the counter slot)
     * @return whether at least one k-mer of the read matched the database
     */
    protected boolean matchRead(final MatcherReadEntry entry, final int index) {
        return matchRead(entry, (MatcherConsumer) getConsumer(index));
    }

    /**
     * Matches one read against the database on behalf of the given consumer, using only that
     * consumer's private state.
     *
     * @param entry the read to match
     * @param consumer the consumer performing the match
     * @return whether the read matched the database at all
     */
    protected boolean matchRead(final MatcherReadEntry entry, final MatcherConsumer consumer) {
        boolean found = false;
        int prints = 0;
        int readTaxErrorCount = taxTree == null ? -1 : 0;

        SmallTaxIdNode taxIdNode;
        int max = entry.readSize - k + 1;
        // Loop-invariant per read: hoisted out of the per-error-k-mer threshold check in the loop.
        double maxReadTaxErrorCountTimesMax = maxReadTaxErrorCount * max;
        SmallTaxIdNode lastTaxid = null;
        int contigLen = 0;
        CountsPerTaxid stats = null;
        // The consumer index is constant for this call, so hoist its per-store-index row once.
        final long[] readNoRow = consumer.readNoPerStat;

        // Resolve all of this read's k-mers up front when the store supports batched lookups, so the
        // loop below consumes results instead of stalling on one cache miss after another.
        final SmallTaxIdNode[] prefetched;
        final long[] prefetchedPositions;
        if (batchStore != null) {
            prefetchKMers(entry, max, consumer);
            prefetched = consumer.prefetchedNodes;
            prefetchedPositions = consumer.prefetchedPositions;
        } else {
            prefetched = null;
            prefetchedPositions = null;
        }

        long kmer = -1;
        long reverseKmer = -1;
        int oldIndex = 0;
        for (int i = 0; i < max; i++) {
            if (kmer == -1) {
                kmer = CGAT.kMerToLongStraight(entry.read, i, k, entry.badPos);
                if (kmer == -1) {
                    oldIndex = i;
                    i = entry.badPos[0];
                } else {
                    reverseKmer = CGAT.kMerToLongReverse(entry.read, i, k, null);
                }
            } else {
                final byte lastBase = entry.read[i + k - 1];
                kmer = CGAT.nextKMerStraight(kmer, lastBase, k);
                if (kmer == -1) {
                    oldIndex = i;
                    i += k - 1;
                } else {
                    reverseKmer = CGAT.nextKMerReverse(reverseKmer, lastBase, k);
                }
            }
            if (kmer == -1) {
                taxIdNode = INVALID_NODE;
            } else if (prefetched != null) {
                taxIdNode = prefetched[i];
                if (taxIdNode != null) {
                    entry.indexPos[0] = prefetchedPositions[i];
                }
            } else {
                taxIdNode = kmerStore.getLong(CGAT.standardKMer(kmer, reverseKmer), entry.indexPos);
                if (taxIdNode != null && uniqueKmerBits != null) {
                    uniqueKmerBits.set(entry.indexPos[0]);
                }
            }
            // Whether this k-mer starts a new contig (its tax node differs from the previous k-mer's).
            // Computed before lastTaxid is updated further below, and used to run the per-contig-only
            // work (the tax-path merge and the stats/reads1KMer resolution) once per contig rather than
            // per k-mer.
            final boolean newContig = taxIdNode != lastTaxid;
            if (readTaxErrorCount != -1) {
                if (taxIdNode == null || taxIdNode == INVALID_NODE) {
                    readTaxErrorCount++;
                    if (maxReadTaxErrorCount >= 0) {
                        if ((maxReadTaxErrorCount >= 1 && readTaxErrorCount > maxReadTaxErrorCount)
                                || (readTaxErrorCount > maxReadTaxErrorCountTimesMax)) {
                            readTaxErrorCount = -1;
                        }
                    }
                } else {
                    // incCount is the per-k-mer vote weight; the tax-path merge is idempotent within a
                    // contig (repeated calls with the same node do not change the path set), so it only
                    // needs to run at the contig start.
                    incCount(taxIdNode, consumer, entry.readNo);
                    if (newContig) {
                        mergeReadTaxidPath(taxIdNode, entry);
                    }
                }
            }
            if (taxIdNode != lastTaxid) {
                if (contigLen > 0) {
                    if (out != null) {
                        printKrakenStyleOut(entry, lastTaxid, contigLen, prints++);
                    }
                    if (stats != null) {
                        synchronized (stats) {
                            // Batched per contig: for a matched contig contigLen equals the number of
                            // its k-mers, so this replaces the former per-k-mer stats.kmers++.
                            stats.kmers += contigLen;
                            stats.contigs++;
                            stats.contigLenSquaredSum += ((long) contigLen) * contigLen;
                            if (contigLen > stats.maxContigLen) {
                                stats.maxContigLen = contigLen;
                                int j = 1;
                                for (; j < entry.readDescriptorSize && j < stats.maxContigDescriptor.length && entry.readDescriptor[j] != ' '; j++) {
                                    stats.maxContigDescriptor[j - 1] = entry.readDescriptor[j];
                                }
                                stats.maxContigDescriptor[j - 1] = 0;
                            }
                        }
                    }
                    contigLen = 0;
                }
            }
            if (taxIdNode == INVALID_NODE) {
                contigLen += i >= max ? max - oldIndex : i - oldIndex + 1;
            }
            else {
                contigLen++;
            }
            lastTaxid = taxIdNode;
            if (taxIdNode != null && taxIdNode != INVALID_NODE) {
                found = true;
                if (newContig) {
                    // 'stats' and the reads1KMer bookkeeping are constant within a contig, so resolve
                    // them once at the contig start; 'stats' is then carried across the contig for the
                    // boundary flush. stats.kmers itself is accumulated per contig in the contig-boundary
                    // block (and the tail), so this hot path no longer locks 'stats' per k-mer.
                    int vi = taxIdNode.getStoreIndex();
                    stats = getCountsPerTaxid(taxIdNode, vi);
                    // reads1KMer is counted once per (read, tax id); the guard row readNoPerCPerStat[index]
                    // is owned by this consumer thread alone, so the check is race-free and only the rare
                    // first hit of a tax id in a read needs the lock.
                    if (readNoRow[vi] != entry.readNo) {
                        readNoRow[vi] = entry.readNo;
                        synchronized (stats) {
                            stats.reads1KMer++;
                        }
                    }
                }
            } else {
                stats = null;
            }
        }
        if (contigLen > 0 && out != null) {
            printKrakenStyleOut(entry, lastTaxid, contigLen, prints);
        }
        if (found) {
            if (contigLen > 0) {
                if (stats != null) {
                    synchronized (stats) {
                        // Batched per contig (final contig): see the boundary block above.
                        stats.kmers += contigLen;
                        stats.contigs++;
                        stats.contigLenSquaredSum += ((long) contigLen) * contigLen;
                        if (contigLen > stats.maxContigLen) {
                            stats.maxContigLen = contigLen;
                            int j = 1;
                            for (; j < entry.readDescriptorSize && j < stats.maxContigDescriptor.length && entry.readDescriptor[j] != ' '; j++) {
                                stats.maxContigDescriptor[j - 1] = entry.readDescriptor[j];
                            }
                            stats.maxContigDescriptor[j - 1] = 0;
                        }
                    }
                }
            }
            if (readTaxErrorCount != -1) {
                int ties = 0;
                for (int i = 0; i < entry.usedPaths; i++) {
                    int sum = sumCounts(entry.readTaxIdNode[i], consumer, entry.readNo);
                    if (sum > entry.counts[0]) {
                        entry.counts[0] = sum;
                        entry.readTaxIdNode[0] = entry.readTaxIdNode[i];
                        ties = 0;
                    } else if (sum == entry.counts[0]) {
                        ties++;
                        entry.counts[ties] = sum;
                        entry.readTaxIdNode[ties] = entry.readTaxIdNode[i];
                    }
                }
                if (threshold > 1) {
                    for (int i = 0; i <= ties; i++) {
                        entry.readTaxIdNode[i] = lowestNodeWhereSumAboveThreshold(entry.readTaxIdNode[i], consumer, entry.readNo, threshold);
                    }
                }
                SmallTaxIdNode node = entry.readTaxIdNode[0];
                for (int i = 1; i <= ties; i++) {
                    node = taxTree.getLowestCommonAncestor(node, entry.readTaxIdNode[i]);
                }
                entry.classNode = node;
                if (node == null) {
                    return false;
                }
                // For 'readKmers', I decided to count in the k-mers from 'entry.readTaxIdNode[0]' and not just 'node'.
                // (They only differ in case of a tie anyways.) But if there is tie, then the k-mers from one of the tie's nodes
                // solidify the LCA in a sense - so the counts from one of the involved paths are included.
                // When threshold > 1, readTaxIdNode[0] was promoted to an ancestor above, so the
                // voting-time entry.counts[0] is stale; recompute sumCounts for the actual node.
                int readKmers = (ties > 0 || threshold > 1)
                        ? sumCounts(entry.readTaxIdNode[0], consumer, entry.readNo) : entry.counts[0];
                int classErrC = max - readKmers;
                if (maxReadClassErrorCount < 0 || (maxReadClassErrorCount >= 1 && classErrC <= maxReadClassErrorCount)
                        || (classErrC <= maxReadClassErrorCount * max)) {
                    double err = ((double) readTaxErrorCount) / max;
                    double classErr = ((double) classErrC) / max;
                    entry.classNode = node;
                    int vi = node.getStoreIndex();
                    if (vi >= 0) {
                        stats = getCountsPerTaxid(node, vi);
                        synchronized (stats) {
                            stats.reads++;
                            stats.readsKmers += readKmers;
                            stats.readsBPs += entry.readSize;
                            stats.errorSum += err;
                            stats.errorSquaredSum += err * err;
                            stats.classErrorSum += classErr;
                            stats.classErrorSquaredSum += classErr * classErr;
                        }
                    }
                    else if (getLogger().isWarnEnabled()) {
                        getLogger().warn("Missing database entry for tax node: " + node);
                    }
                }
            }
        }

        return found;
    }

    /**
     * Returns the statistics object for the given store index {@code vi}, lazily creating
     * it for the given node in a thread-safe way if necessary.
     *
     * @param node the tax id node the statistics belong to
     * @param vi   the store index selecting the statistics slot
     * @return the (possibly newly created) statistics object for the node
     */
    protected final CountsPerTaxid getCountsPerTaxid(final SmallTaxIdNode node, final int vi) {
        CountsPerTaxid stats = statsIndex[vi];
        if (stats == null) {
            synchronized (statsIndex) {
                if (statsIndex[vi] == null) {
                    statsIndex[vi] = new CountsPerTaxid(node.getLevel(), node.getTaxId(), initialReadSize);
                }
                stats = statsIndex[vi];
            }
        }
        return stats;
    }

    /**
     * Merges the given tax id node into the read's set of candidate taxonomic paths, collapsing paths
     * that are ancestors of one another so that only the most specific nodes are kept. This operation
     * is idempotent for a node already represented in the path set, so the matcher runs it only once
     * per contig rather than once per k-mer.
     *
     * @param node  the tax id node hit by a k-mer of the current read
     * @param entry the read together with its per-read working state
     */
    // Made final for potential inlining by JVM
    protected final void mergeReadTaxidPath(final SmallTaxIdNode node, final MatcherReadEntry entry) {
        boolean found = false;
        for (int i = 0; i < entry.usedPaths; i++) {
            if (taxTree.isAncestorOf(node, entry.readTaxIdNode[i])) {
                entry.readTaxIdNode[i] = node;
                found = true;
                break;
            } else if (taxTree.isAncestorOf(entry.readTaxIdNode[i], node)) {
                found = true;
                break;
            }
        }
        if (!found) {
            if (entry.usedPaths < maxPaths) {
                entry.readTaxIdNode[entry.usedPaths] = node;
                entry.usedPaths++;
            }
        }
    }

    /**
     * Appends a single Kraken-style segment ({@code taxid:contigLen}) to the read's output
     * buffer, using {@code A} for invalid (ambiguous) and {@code 0} for unmatched k-mers.
     *
     * @param entry     the read whose output buffer is appended to
     * @param taxid     the tax id node of the segment, or {@code null}/{@code INVALID_NODE}
     * @param contigLen the length of the contiguous k-mer run
     * @param state     the segment index (a leading space is added when non-zero)
     */
    protected void printKrakenStyleOut(final MatcherReadEntry entry, final SmallTaxIdNode taxid, final int contigLen, final int state) {
        if (state != 0) {
            entry.printChar(' ');
        }
        if (taxid == INVALID_NODE) {
            entry.printChar('A');
        }
        else if (taxid == null) {
            entry.printChar('0');
        } else {
            entry.printString(taxid.getTaxId());
        }
        entry.printChar(':');
        entry.printInt(contigLen);
    }

    /**
     * Per-read working state used during matching: the candidate taxonomic paths and
     * their counts, the resolved classification node, and a growable byte buffer holding
     * the read's Kraken-style output.
     */
    public static class MatcherReadEntry extends ReadEntry {
        private final static byte[] U = new byte[] { 'U', '\t'};
        private final static byte[] C = new byte[] { 'C', '\t'};
        private final static int CLASS_TAX_BUFFER_SIZE = 128;

        /** Growable byte buffer holding the read's Kraken-style output. */
        public byte[] buffer;
        /** Current write position within {@link #buffer}. */
        public int bufferPos;

        private byte[] classTaxBuffer;

        /** Single-element scratch array receiving the position of a bad (non-CGAT) base. */
        public int[] badPos = new int[1];

        /** Number of candidate taxonomic paths currently used. */
        public int usedPaths;
        /** The candidate taxonomic path nodes tracked for this read. */
        public SmallTaxIdNode[] readTaxIdNode;
        /** Vote counts associated with the candidate paths. */
        public int[] counts;
        /** Single-element scratch array holding the store index of the last matched k-mer. */
        public long[] indexPos;
        /** The resolved classification node for this read, or {@code null} if unclassified. */
        public SmallTaxIdNode classNode;

        /**
         * Creates a read entry with buffers sized for the given number of candidate paths.
         *
         * @param maxReadSizeBytes the initial read size in bytes
         * @param withProbs        whether quality/probability information is processed
         * @param paths            the maximum number of candidate taxonomic paths
         */
        public MatcherReadEntry(int maxReadSizeBytes, boolean withProbs, int paths) {
            super(maxReadSizeBytes, withProbs);

            buffer = null;
            readTaxIdNode = new SmallTaxIdNode[paths];
            counts = new int[paths];
            indexPos = new long[1];
        }

        /**
         * Appends a single character (as a byte) to the output buffer.
         *
         * @param c the character to append
         */
        public void printChar(final char c) {
            growPrintBuffer(1);
            buffer[bufferPos++] = (byte) c;
        }

        /**
         * Appends the bytes of the given string to the output buffer.
         *
         * @param s the string to append
         */
        public void printString(final String s) {
            int len = s.length();
            growPrintBuffer(len);
            s.getBytes(0, len, buffer, bufferPos);
            bufferPos += len;
        }

        /**
         * Ensures the output buffer has room for at least {@code additionalSize} more
         * bytes, allocating or doubling it as needed.
         *
         * @param additionalSize the number of additional bytes that must fit
         */
        protected void growPrintBuffer(int additionalSize) {
            if (buffer == null) {
                buffer = new byte[additionalSize + 1024];
                return;
            }

            int newLen = buffer.length;
            while (bufferPos + additionalSize > newLen) {
                newLen *= 2;
            }
            if (newLen > buffer.length) {
                byte[] newBuffer = new byte[newLen];
                System.arraycopy(buffer, 0, newBuffer, 0, bufferPos);
                buffer = newBuffer;
            }
        }

        /**
         * Appends the decimal representation of the given integer to the output buffer.
         *
         * @param value the integer value to append
         */
        public void printInt(final int value) {
            growPrintBuffer(11); // Decimal version of int has 11 bytes max.
            bufferPos = ByteArrayUtil.intToByteArray(value, buffer, bufferPos);
        }

        /**
         * Writes this read's Kraken-style classification line to the given stream: the
         * classified/unclassified flag, read descriptor, assigned tax id, read length and
         * the buffered per-segment output.
         *
         * @param out the stream to write the classification line to
         * @throws IOException if writing to the stream fails
         */
        public void writeMatchDetails(OutputStream out) throws IOException {
            if (buffer == null) {
                return;
            }
            if (classNode == null) {
                out.write(U);
            } else {
                out.write(C);
            }
            int index = ByteArrayUtil.indexOf(readDescriptor, 1, readDescriptorSize, ' ');
            // ByteArrayUtil.println(readDescriptor, 1, index == -1 ? readDescriptorSize : index, System.out);
            out.write(readDescriptor, 1, (index == -1 ? readDescriptorSize : index) - 1);
            out.write('\t');
            if (classNode == null) {
                out.write('0');
            } else {
                if (classTaxBuffer == null) {
                    classTaxBuffer = new byte[CLASS_TAX_BUFFER_SIZE];
                }
                String tax = classNode.getTaxId();
                int len =  tax.length();
                tax.getBytes(0, len, classTaxBuffer, 0);
                out.write(classTaxBuffer, 0, len);
            }
            out.write('\t');
            if (classTaxBuffer == null) {
                classTaxBuffer = new byte[CLASS_TAX_BUFFER_SIZE];
            }
            int len =  ByteArrayUtil.intToByteArray(readSize, classTaxBuffer, 0);
            out.write(classTaxBuffer, 0, len);
            out.write('\t');
            out.write(buffer, 0, bufferPos);
            out.write('\n');
        }
   }

   /**
    * Callback invoked after each read has been matched, receiving the read entry and
    * whether it matched.
    */
   public interface AfterMatchCallback {
        /**
         * Invoked after a read has been matched.
         *
         * @param entry the matched read entry
         * @param found whether the read matched at least one k-mer
         */
        void afterMatch(MatcherReadEntry entry, boolean found);
   }
}
