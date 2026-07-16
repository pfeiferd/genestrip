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
import org.metagene.genestrip.fastq.AbstractLoggingFastqStreamer;
import org.metagene.genestrip.io.StreamProvider;
import org.metagene.genestrip.io.StreamingResource;
import org.metagene.genestrip.io.StreamingResourceListStream;
import org.metagene.genestrip.io.StreamingResourceStream;
import org.metagene.genestrip.store.KMerStore;
import org.metagene.genestrip.store.KMerUniqueCounterBits;
import org.metagene.genestrip.tax.SmallTaxTree;
import org.metagene.genestrip.tax.SmallTaxTree.SmallTaxIdNode;
import org.metagene.genestrip.util.ByteArrayUtil;
import org.metagene.genestrip.util.CGAT;

import it.unimi.dsi.fastutil.objects.Object2LongMap;

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
    /** Maximum number of top k-mer counts to retain per tax id. */
    protected final int maxKmerResCounts;

    // Turned from KMerUniqueCounter to KMerUniqueCounterBits for potential method inlining.
    /** Optional counter of unique k-mers per tax id, or {@code null} if unique counting is disabled. */
    protected KMerUniqueCounterBits uniqueCounter;
    /** Per-store-index statistics accumulators, indexed by a tax id node's store index. */
    protected final CountsPerTaxid[] statsIndex;
    /** Per-consumer, per-store-index last read number seen, used to count reads with at least one k-mer. */
    protected final long readNoPerCPerStat[][];

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
     * @param maxKmerResCounts      the maximum number of top k-mer counts to retain per tax id
     * @param taxTree               the taxonomy tree used for voting and LCA resolution
     * @param maxPaths              the maximum number of candidate taxonomic paths per read
     * @param maxReadTaxErrorCount  the maximum allowed tax error count per read
     * @param maxReadClassErrorCount the maximum allowed classification error count per read
     * @param writeAll              whether to write all reads, not just classified ones
     * @param threshold             the minimum vote count threshold for classification
     * @param dbMD5                 the MD5 checksum identifying the database
     */
    public FastqKMerMatcher(KMerStore<SmallTaxIdNode> kmerStore, int initialReadSize, int maxQueueSize,
                            ExecutionContext bundle, boolean withProbs, int maxKmerResCounts, SmallTaxTree taxTree, int maxPaths,
                            double maxReadTaxErrorCount, double maxReadClassErrorCount, boolean writeAll, int threshold, String dbMD5) {
        super(kmerStore.getK(), initialReadSize, maxQueueSize, bundle, withProbs, maxPaths);
        consumers = bundle.getThreads() <= 0 ? 1 : bundle.getThreads();
        this.kmerStore = kmerStore;
        this.statsIndex = new CountsPerTaxid[kmerStore.getNValues()];
        this.readNoPerCPerStat = new long[consumers][kmerStore.getNValues()];
        this.initialReadSize = initialReadSize;
        this.maxKmerResCounts = maxKmerResCounts;
        this.taxTree = taxTree;
        this.maxReadTaxErrorCount = maxReadTaxErrorCount;
        this.maxReadClassErrorCount = maxReadClassErrorCount;
        this.maxPaths = maxPaths;
        this.writeAll = writeAll;
        this.threshold = threshold;
        this.dbMD5 = dbMD5;
        if (taxTree != null) {
            taxTree.initCountSize(consumers);
        }
    }

    @Override
    protected ReadEntry createReadEntry(int initialReadSizeBytes, boolean withProbs, Object... config) {
        return new MatcherReadEntry(initialReadSizeBytes, withProbs, (int) config[0]);
    }

    /**
     * Convenience overload of {@link #runMatcher(StreamingResourceStream, File, File, KMerUniqueCounterBits)}
     * that matches a single FASTQ resource.
     *
     * @param fastq             the FASTQ resource to process
     * @param filteredFile      optional file to which matched reads are written, or {@code null}
     * @param krakenOutStyleFile optional file for Kraken-style per-read output, or {@code null}
     * @param uniqueCounter     optional unique-k-mer counter, or {@code null} to skip unique counting
     * @return the aggregated matching result
     * @throws IOException if reading the FASTQ resource or writing the output files fails
     */
    public MatchingResult runMatcher(StreamingResource fastq, File filteredFile, File krakenOutStyleFile,
                                     KMerUniqueCounterBits uniqueCounter) throws IOException {
        return runMatcher(new StreamingResourceListStream(fastq), filteredFile, krakenOutStyleFile, uniqueCounter);
    }

    /**
     * Matches all reads of the given FASTQ streams against the database and returns the
     * aggregated result.
     *
     * @param fastqs            the FASTQ resources to process
     * @param filteredFile      optional file to which matched reads are written, or {@code null}
     * @param krakenOutStyleFile optional file for Kraken-style per-read output, or {@code null}
     * @param uniqueCounter     optional unique-k-mer counter, or {@code null} to skip unique counting
     * @return the aggregated matching result
     * @throws IOException if reading the FASTQ streams or writing the output files fails
     */
    public MatchingResult runMatcher(StreamingResourceStream fastqs, File filteredFile, File krakenOutStyleFile,
                                     KMerUniqueCounterBits uniqueCounter) throws IOException {
        try (OutputStream lindexed = filteredFile != null ? StreamProvider.getOutputStreamForFile(filteredFile) : null;
             // A PrintStream is implicitly synchronized. So we don't need to worry about
             // multi threading when using it.
             PrintStream lout = krakenOutStyleFile != null
                     ? new PrintStream(StreamProvider.getOutputStreamForFile(krakenOutStyleFile), false, StandardCharsets.UTF_8)
                     : null) {
            indexed = lindexed;
            out = lout;

            initStats();
            initUniqueCounter(uniqueCounter);
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

        Map<String, short[]> countMap = null;
        if (uniqueCounter != null) {
            Object2LongMap<String> counts = uniqueCounter.getUniqueKmerCounts();
            for (CountsPerTaxid stats : statsIndex) {
                if (stats != null) {
                    stats.uniqueKmers = counts.getLong(stats.getTaxid());
                }
            }
            if (uniqueCounter instanceof KMerUniqueCounterBits) {
                if (((KMerUniqueCounterBits) uniqueCounter).isWithCounts()) {
                    countMap = ((KMerUniqueCounterBits) uniqueCounter).getMaxCountsCounts(maxKmerResCounts);
                    for (CountsPerTaxid stats : statsIndex) {
                        if (stats != null) {
                            stats.maxKMerCounts = countMap.get(stats.getTaxid());
                        }
                    }
                }
            }
            this.uniqueCounter = null;
        } else {
            for (CountsPerTaxid stats : statsIndex) {
                if (stats != null) {
                    stats.uniqueKmers = -1;
                }
            }
        }

        return new MatchingResult(kmerStore.getK(), taxid2Stats, dbMD5, totalReads, totalKMers, totalBPs,
                countMap == null ? null : countMap.get(null));
    }

    @Override
    protected void readFastq(InputStream inputStream, boolean fasta) throws IOException {
        try {
            if (taxTree != null) {
                taxTree.resetCounts(this);
            }
            for (long[] a : readNoPerCPerStat) {
                Arrays.fill(a, -1);
            }
            super.readFastq(inputStream, fasta);
        } finally {
            if (taxTree != null) {
                taxTree.releaseOwner();
            }
        }
    }

    // Package private for testing purposes.
    void initStats() {
        Arrays.fill(statsIndex, null);
    }

    void initUniqueCounter(KMerUniqueCounterBits uniqueCounter) {
        // Turned from KMerUniqueCounter to KMerUniqueCounterBits for potential method inlining.
        this.uniqueCounter = uniqueCounter;
        if (uniqueCounter != null) {
            uniqueCounter.clear();
        }
    }

    @Override
    // Made final for potential inlining by JVM
    protected final void nextEntry(ReadEntry entry, int index) throws IOException {
        MatcherReadEntry myEntry = (MatcherReadEntry) entry;
        myEntry.bufferPos = 0;

        myEntry.usedPaths = 0;
        myEntry.classNode = null;
        for (int i = 0; i < maxPaths; i++) {
            myEntry.readTaxIdNode[i] = null;
            myEntry.counts[i] = 0;
        }

        boolean found = matchRead(myEntry, index);
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
        final long[] readNoRow = readNoPerCPerStat[index];

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
            taxIdNode = kmer == -1 ? INVALID_NODE :
                    kmerStore.getLong(CGAT.standardKMer(kmer, reverseKmer), entry.indexPos);
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
                    taxTree.incCount(taxIdNode, index, entry.readNo);
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
                if (uniqueCounter != null) {
                    // This is a considerable optimization as found via profiling:
                    // Old version:
                    // uniqueCounter.put(CGAT.standardKMer(kmer, reverseKmer), taxIdNode.getTaxId(), entry.indexPos[0]);
                    // Faster version:
                    uniqueCounter.putInlined(entry.indexPos[0]);
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
                    int sum = taxTree.sumCounts(entry.readTaxIdNode[i], index, entry.readNo);
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
                        entry.readTaxIdNode[i] = taxTree.lowestNodeWhereSumAboveThreshold(entry.readTaxIdNode[i], index, entry.readNo, threshold);
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
                        ? taxTree.sumCounts(entry.readTaxIdNode[0], index, entry.readNo) : entry.counts[0];
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
