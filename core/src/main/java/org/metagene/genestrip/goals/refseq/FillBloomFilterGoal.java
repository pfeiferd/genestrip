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
package org.metagene.genestrip.goals.refseq;

import java.io.File;
import java.io.IOException;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.metagene.genestrip.ExecutionContext;
import org.metagene.genestrip.GSConfigKey;
import org.metagene.genestrip.GSGoalKey;
import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.probfilter.BloomFilter;
import org.metagene.genestrip.probfilter.ProbFilter;
import org.metagene.genestrip.probfilter.XORBloomFilter;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.refseq.AbstractRefSeqFastaReader;
import org.metagene.genestrip.refseq.AbstractStoreFastaReader;
import org.metagene.genestrip.refseq.AccessionMap;
import org.metagene.genestrip.refseq.RefSeqCategory;
import org.metagene.genestrip.refseq.ReworkingStoreFastaReader;
import org.metagene.genestrip.store.RadixKMerStore;
import org.metagene.genestrip.tax.Rank;
import org.metagene.genestrip.tax.TaxTree;
import org.metagene.genestrip.tax.TaxTree.IDStringGenerator;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;
import org.metagene.genestrip.util.ByteArrayUtil;
import org.metagene.genestrip.util.StringLongDigitTrie;
import org.metagene.genestrip.tax.TaxNodeSelection;

/**
 * Goal that estimates the deduplicated database size: it streams all selected k-mers through a
 * temporary Bloom filter and counts the distinct k-mers, in total and per {@link RadixKMerStore}
 * radix bucket, producing a {@link DBSize} used to size the final k-mer store.
 *
 * @param <P> the project type
 */
public class FillBloomFilterGoal<P extends GSProject> extends FastaReaderGoal<FillBloomFilterGoal.DBSize, P> implements Goal.LogHeapInfo {
    /**
     * Combined deduplicated DB size: the total estimated number of distinct k-mers and, when
     * available, the per-{@link RadixKMerStore} radix-bucket counts (length {@code 2^radixStoreBits},
     * suitable as the {@code bucketSizes} argument of a {@link RadixKMerStore} created with the same
     * {@code radixStoreBits}). {@code bucketSizes} is {@code null} for sizing strategies that do not
     * compute per-bucket counts.
     */
    public static class DBSize implements Serializable {
        private static final long serialVersionUID = 1L;

        /** The per-radix-bucket distinct k-mer counts, or {@code null} if not computed. */
        private final int[] bucketSizes;
        /**
         * The set of distinct store values (tax id strings) that the fill will use, collected during
         * the counting pass so they can be registered up front, or {@code null} if not collected. Only
         * complete when the fill does not synthesize artificial data/file/id nodes (those are created
         * during the fill itself and cannot be pre-collected).
         */
        private final Set<String> values;
        /**
         * Creates a DB size holder.
         *
         * @param bucketSizes the per-radix-bucket counts, or {@code null} if not computed
         * @param values the distinct fill values collected up front, or {@code null} if not collected
         */
        public DBSize(int[] bucketSizes, Set<String> values) {
            this.bucketSizes = bucketSizes;
            this.values = values;
        }

        /**
         * Returns the total estimated number of distinct k-mers.
         *
         * @return the total estimated number of distinct k-mers
         */
        public long getSize() {
            long sum = 0;
            for (int i = 0;  i < bucketSizes.length; i++) {
                sum += bucketSizes[i];
            }
            return sum;
        }

        /**
         * Returns the per-radix-bucket distinct k-mer counts.
         *
         * @return the per-radix-bucket counts, or {@code null} if not computed
         */
        public int[] getBucketSizes() {
            return bucketSizes;
        }

        /**
         * Returns the distinct fill values collected up front (for pre-registering the store's value
         * map), or {@code null} if they were not collected.
         *
         * @return the distinct fill values, or {@code null}
         */
        public Set<String> getValues() {
            return values;
        }
    }

    private final ObjectGoal<AccessionMap, P> accessionMapGoal;
    private final ObjectGoal<TaxTree, P> taxTreeGoal;
    private final ObjectGoal<FillSizeGoal.KMerCounts, P> sizeGoal;

    // The temporary size-estimation filter is always an XOR/Murmur filter. The store sizing is
    // derived from the readers' per-radix-bucket counts of distinct k-mers, not from the filter.
    /**
     * The false-positive rate beyond which the k-mer counts are refused rather than used. Set from
     * where {@link BloomFilter#estimateDistinctValues(long)} was measured to leave the range in which
     * it is good to a fraction of a percent: at a reached rate of 0.14 it is out by 0.3 %, at 0.33 by
     * 3 %, at 0.46 by 9 %. Refusing above 0.25 therefore draws the line where the number stops being
     * an estimate, rather than where it becomes nonsense.
     */
    private static final double MAX_TRUSTED_FPP = 0.25;


    private ProbFilter filter;
    // Number of low k-mer bits used as the radix (from config); the radix store created later must
    // use the same value. Set in doMakeThis().
    private int radixBits;
    // Each reader thread counts its per-radix-bucket k-mers and collects its distinct store values into
    // its own thread-local arrays/sets (no synchronization); doMakeThis() merges them after readFastas().
    private final List<MyFastaReader> readers = new ArrayList<>();

    /**
     * Creates the goal, wiring the accession-map and expected-size goals alongside the FASTA inputs.
     *
     * @param project the project type
     * @param bundle the execution context providing threading and shared services
     * @param categoriesGoal the goal supplying the selected RefSeq categories
     * @param taxNodesGoal the goal supplying the selected taxonomic nodes
     * @param fnaFilesGoal the goal supplying the downloaded RefSeq FASTA files
     * @param additionalGoal the goal supplying additional FASTA files mapped to tax nodes
     * @param accessionMapGoal the goal supplying the accession-to-tax-id map
     * @param taxTreeGoal the goal supplying the taxonomy tree (into which artificial fill nodes are created)
     * @param sizeGoal the goal supplying the k-mer counts the filter is sized from
     * @param deps the additional goals this goal depends on
     */
    @SafeVarargs
    public FillBloomFilterGoal(P project, ExecutionContext bundle, ObjectGoal<Set<RefSeqCategory>, P> categoriesGoal,
                               ObjectGoal<TaxNodeSelection, P> taxNodesGoal, RefSeqFnaFilesDownloadGoal fnaFilesGoal,
                               ObjectGoal<Map<File, TaxIdNode>, P> additionalGoal,
                               ObjectGoal<AccessionMap, P> accessionMapGoal, ObjectGoal<TaxTree, P> taxTreeGoal,
                               ObjectGoal<FillSizeGoal.KMerCounts, P> sizeGoal, Goal<P>... deps) {
        super(project, GSGoalKey.TEMPINDEX, bundle, categoriesGoal, taxNodesGoal, fnaFilesGoal, additionalGoal, Goal.append(deps, accessionMapGoal, taxTreeGoal, sizeGoal));
        this.accessionMapGoal = accessionMapGoal;
        this.taxTreeGoal = taxTreeGoal;
        this.sizeGoal = sizeGoal;
    }

    /**
     * Eagerly cleans this goal's result once all dependent goals have been made, to free memory.
     */
    @Override
    protected void allDependentsMade() {
        // To save memory...
        doCleanThis();
    }

    @Override
    protected void doMakeThis() {
        FillSizeGoal.KMerCounts counts = sizeGoal.get();
        GSConfigKey.BloomFilterSizing sizing = (GSConfigKey.BloomFilterSizing) configValue(GSConfigKey.BLOOM_FILTER_SIZING);
        // AUTO is DISTINCT with a second attempt held in reserve. The reserve is used at most once:
        // the count with duplicates is exact and an upper bound, so a filter built for it cannot be
        // too small and there is nothing a third attempt could try.
        boolean fromDistinct = sizing != GSConfigKey.BloomFilterSizing.UPPER_BOUND;
        boolean mayRetry = sizing == GSConfigKey.BloomFilterSizing.AUTO;
        try {
            try {
                set(onePass(counts, fromDistinct));
            } catch (UnusableFilterException e) {
                if (!mayRetry) {
                    throw new IllegalStateException(e.explain(counts, fromDistinct,
                            GSConfigKey.BloomFilterSizing.UPPER_BOUND.getName()));
                }
                if (getLogger().isWarnEnabled()) {
                    getLogger().warn(e.explain(counts, fromDistinct, null)
                            + " Reading the sequences again with the exact count, as '"
                            + GSConfigKey.BLOOM_FILTER_SIZING.getName() + "="
                            + GSConfigKey.BloomFilterSizing.AUTO.getName() + "' asks for. This doubles the"
                            + " time this goal takes; set it to '"
                            + GSConfigKey.BloomFilterSizing.UPPER_BOUND.getName()
                            + "' to go straight there next time.");
                }
                // The first pass has ended its own consumers, so only this goal's own state is
                // left to reset before another one reads.
                readers.clear();
                try {
                    set(onePass(counts, false));
                } catch (UnusableFilterException second) {
                    // The exact count is an upper bound, so this cannot happen for want of size, and
                    // reporting it as though it could would send the reader looking in the wrong place.
                    throw new IllegalStateException(second.explain(counts, false, null)
                            + " This second attempt was sized from the exact count including duplicates,"
                            + " which no data can exceed, so the cause lies elsewhere - the k-mers counted"
                            + " here outnumber every k-mer the size goal reported seeing.");
                }
            }
        } catch (IOException e) {
            throw new RuntimeException(e);
        } finally {
            filter = null;
            readers.clear();
            // The filter is sized for the whole database and is dead the moment this goal is done
            // with it, while the goal that follows is the one that allocates the k-mer store. Handing
            // the memory back here rather than waiting for the virtual machine to feel the pressure is
            // therefore worth the collection - this is about the memory, not about the figure the
            // logging takes, which collects on its own account and only when it is switched on.
            System.gc();
        }
    }

    /**
     * Reads all sequences once through a temporary filter of the given sizing and returns the
     * database size derived from what it counted.
     *
     * @param counts the two k-mer counts the filter may be sized from
     * @param fromDistinct whether to size it from the estimated distinct count rather than the exact one
     * @return the resulting database size
     * @throws UnusableFilterException if the filter ended up too full for its counts to be trusted
     * @throws IOException if a sequence file cannot be read
     */
    protected DBSize onePass(FillSizeGoal.KMerCounts counts, boolean fromDistinct)
            throws UnusableFilterException, IOException {
        double tempFpp = doubleConfigValue(GSConfigKey.TEMP_BLOOM_FILTER_FPP);
        long sizingCount = fromDistinct ? counts.getDistinct() : counts.getWithDuplicates();
        if (getLogger().isInfoEnabled()) {
            getLogger().info("Sizing the temporary bloom filter for " + sizingCount + " k-mers, the "
                    + (fromDistinct ? "estimated distinct" : "exact with-duplicates") + " count of " + counts);
        }
        // The temporary size-estimation filter is filled concurrently by the reader threads, so it
        // must support a thread-safe putLong; the two used here do (via their bit vector's bucket locks),
        // as does BlockedBloomFilter, which is nonetheless not offered as an option here.
        // Kept as a BloomFilter beside the field, which is a ProbFilter: the correction below
        // asks it for the false-positive rate it reached, which only a bloom filter can answer.
        BloomFilter tempFilter = booleanConfigValue(GSConfigKey.XOR_BLOOM_HASH) ?
                new XORBloomFilter(tempFpp, sizingCount) :
                new BloomFilter(tempFpp, sizingCount);
        filter = tempFilter;
        radixBits = intConfigValue(GSConfigKey.RADIX_STORE_BITS);
        logHeapInfo();
        readFastas();
        // Merge the per-thread counts/values now that all readers have finished (no synchronization
        // needed during the read).
        int[] bucketSizes = new int[1 << radixBits];
        Set<String> collectedValues = new HashSet<>();
        for (MyFastaReader reader : readers) {
            int[] readerBuckets = reader.getBucketSizes();
            for (int i = 0; i < bucketSizes.length; i++) {
                bucketSizes[i] += readerBuckets[i];
            }
            collectedValues.addAll(reader.getValues());
        }
        // We have to account for the missing entries in the bloom filter due to inherent FPP: a k-mer
        // the filter falsely reports as present is never counted, so the counts are short by about the
        // filter's false-positive rate. Correcting for it lets 'bloomFilterGoal' run at a low FPP and
        // saves memory during db construction.
        //
        // The rate to correct by is neither the one the filter was sized for nor the one it ends up
        // at. A k-mer is offered while the filter is still filling, so what hides it is the rate at
        // that moment: with c counted so far, dc/dT = 1 - fpp(c), and the correction is the integral
        // of 1/(1 - fpp(c)) over c, which BloomFilter.estimateDistinctValues does.
        long entries = 0;
        for (int i = 0; i < bucketSizes.length; i++) {
            entries += bucketSizes[i];
        }
        double reachedFpp = tempFilter.getFpp(entries);
        if (reachedFpp > MAX_TRUSTED_FPP) {
            // Beyond this rate the filter hides k-mers faster than the correction recovers them, so
            // the counts are a lower bound of unknown tightness. Sizing the k-mer store from them
            // would not make a smaller database but a wrong one.
            throw new UnusableFilterException(sizingCount, entries, reachedFpp);
        }
        double correctedEntries = tempFilter.estimateDistinctValues(entries);
        // Applied as one factor across the buckets, so that they keep their proportions and their sum
        // stays the corrected total.
        double factor = entries > 0 ? correctedEntries / entries : 1d;
        int[] correctedBucketSizes = new int[bucketSizes.length];
        for (int i = 0; i < correctedBucketSizes.length; i++) {
            correctedBucketSizes[i] = (int) (bucketSizes[i] * factor) + 1;
        }
        // The genome selection travels with the value set: this pass created the artificial nodes those
        // values name, so the fill has to admit exactly the genomes admitted here. It is the sizing
        // pass's selection, already frozen there and handed on unchanged - this pass followed it rather
        // than making one of its own, and the fill will do the same.
        DBSize dbSize = new DBSize(correctedBucketSizes, collectedValues);
        if (getLogger().isInfoEnabled()) {
            getLogger().info("Bloom filter size in kmers: " + dbSize.getSize());
            getLogger().info("Duplication factor: " + ((double) counts.getWithDuplicates()) / dbSize.getSize());
        }
        return dbSize;
    }

    /**
     * Signals that a pass ended with a filter too full for what it counted to be worth using.
     */
    protected static class UnusableFilterException extends Exception {
        private static final long serialVersionUID = 1L;

        private final long sizingCount;
        private final long entries;
        private final double reachedFpp;

        UnusableFilterException(long sizingCount, long entries, double reachedFpp) {
            super("temporary bloom filter too full to count with");
            this.sizingCount = sizingCount;
            this.entries = entries;
            this.reachedFpp = reachedFpp;
        }

        /**
         * Returns the numbers of this failure in words, and what to do about it.
         *
         * @param counts the two counts the filter could have been sized from
         * @param fromDistinct whether the failed pass was sized from the estimate
         * @param remedyValue the configuration value to recommend, or {@code null} for no recommendation
         * @return the explanation
         */
        String explain(FillSizeGoal.KMerCounts counts, boolean fromDistinct, String remedyValue) {
            StringBuilder text = new StringBuilder();
            text.append("The temporary bloom filter is too full for its counts to be trusted: it was sized for ")
                    .append(sizingCount).append(" k-mers")
                    .append(fromDistinct ? " (the estimated number of distinct ones)"
                            : " (the exact number including duplicates)")
                    .append(" but holds ").append(entries)
                    .append(", which brings its false-positive rate to ")
                    .append(String.format("%.3f", reachedFpp))
                    .append(", where anything above ").append(MAX_TRUSTED_FPP)
                    .append(" is unusable. A database size derived from it would be too small, and the k-mer")
                    .append(" store built to that size would lose k-mers.");
            if (fromDistinct) {
                text.append(" The sizing rests on a HyperLogLog estimate, and this one fell short of the ")
                        .append(entries).append(" k-mers actually counted.");
                if (remedyValue != null) {
                    text.append(" Set '").append(GSConfigKey.BLOOM_FILTER_SIZING.getName()).append("=")
                            .append(remedyValue)
                            .append("' in the project's config.properties and run the goal again: the count")
                            .append(" including duplicates is exact and an upper bound, so it cannot fall short")
                            .append(" - at the price of a filter larger by the duplication factor of the data (")
                            .append(String.format("%.1f", counts.getWithDuplicates()
                                    / (double) Math.max(1, counts.getDistinct())))
                            .append(" here, i.e. ").append(counts.getWithDuplicates())
                            .append(" k-mers to size it for).");
                }
            }
            return text.toString();
        }
    }

    @Override
    protected void afterReadFastas(AbstractRefSeqFastaReader.StringLong2DigitTrie contigsPerTaxid) {
        if (getLogger().isDebugEnabled()) {
            List<StringLongDigitTrie.StringLong> list = new ArrayList<>();
            contigsPerTaxid.collect(list);
            getLogger().debug("Contigs ber taxid:");
            getLogger().debug(list);
        }
    }

    @Override
    protected AbstractStoreFastaReader createFastaReader(AbstractRefSeqFastaReader.StringLong2DigitTrie contigsPerTaxid) {
        // Per-reader artificial-tax-id generator (its buffer is mutated), mirroring FillDBGoal. The
        // shared artificial counter lives on the tree, so ids are unique across reader threads.
        byte[] idBuffer = new byte[128];
        idBuffer[0] = '0';
        idBuffer[1] = '0';
        IDStringGenerator idStringGenerator = counter -> {
            int len = ByteArrayUtil.intToByteArray(counter, idBuffer, 2);
            return new String(idBuffer, 0, len);
        };
        MyFastaReader fastaReader = new MyFastaReader(intConfigValue(GSConfigKey.FASTA_LINE_SIZE_BYTES),
                taxNodesGoal.get().getSelected(),
                isIncludeRefSeqFna() ? accessionMapGoal.get() : null,
                intConfigValue(GSConfigKey.KMER_SIZE),
                filter,
                intConfigValue(GSConfigKey.MAX_DUST),
                intConfigValue(GSConfigKey.KMER_SAMPLING),
                booleanConfigValue(GSConfigKey.ASSEMBLY_ACCESSIONS_ONLY),
                contigsPerTaxid,
                booleanConfigValue(GSConfigKey.ENABLE_LOWERCASE_BASES),
                taxTreeGoal.get(),
                booleanConfigValue(GSConfigKey.DATA_NODES),
                booleanConfigValue(GSConfigKey.FILE_NODES),
                booleanConfigValue(GSConfigKey.ID_NODES),
                idStringGenerator,
                (Rank) configValue(GSConfigKey.FOLD_TAXA_BELOW));
        readers.add(fastaReader);
        return fastaReader;
    }

    /**
     * FASTA reader that adds each not-yet-seen k-mer to the temporary Bloom filter and counts the
     * distinct k-mers per radix bucket in a thread-safe way.
     */
    protected class MyFastaReader extends ReworkingStoreFastaReader {
        private final ProbFilter filter;
        // Thread-local per-radix-bucket k-mer counts and distinct store values - this reader runs on a
        // single thread, so no synchronization is needed; doMakeThis() merges all readers afterwards.
        private final int[] bucketSizes = new int[1 << radixBits];
        private final Set<String> values = new HashSet<>();
        // Last node whose taxid was collected, so the (rare) contig-boundary collection skips the
        // repeated per-k-mer adds within a contig.
        private TaxIdNode lastCollectedNode;

        /**
         * Creates the reader that fills the temporary Bloom filter and counts distinct k-mers per
         * radix bucket.
         *
         * @param bufferSize the FASTA line read-buffer size in bytes
         * @param taxNodes the taxonomic nodes to keep k-mers for
         * @param accessionMap the accession-to-tax-id map, or {@code null} if not used
         * @param k the k-mer size
         * @param filter the temporary Bloom filter to fill
         * @param maxDust the maximum allowed low-complexity (dust) run length
         * @param kMerSampling the k-mer sampling step size
         * @param assemblyAccessionsOnly whether only genomic accessions are considered, dropping `NG_`, `NT_` and `NW_`
         * @param contigsPerTaxid the per-tax-id contig counter
         * @param enableLowerCaseBases whether lower-case bases are processed
         * @param taxTree the taxonomy tree into which artificial nodes are created
         * @param dataNodes whether to rework into an artificial {@code DATA} node
         * @param fileNodes whether to rework into an artificial {@code FILE} node
         * @param idNodes whether to rework into an artificial {@code ID} node
         * @param idStringGenerator generator for artificial tax ids
         */
        public MyFastaReader(int bufferSize, Set<TaxIdNode> taxNodes, AccessionMap accessionMap, int k,
                             ProbFilter filter, int maxDust, int kMerSampling, boolean assemblyAccessionsOnly, StringLong2DigitTrie contigsPerTaxid, boolean enableLowerCaseBases,
                             TaxTree taxTree, boolean dataNodes, boolean fileNodes, boolean idNodes, IDStringGenerator idStringGenerator,
                             Rank foldTaxaBelow) {
            super(bufferSize, taxNodes, accessionMap, k, maxDust, kMerSampling, assemblyAccessionsOnly, contigsPerTaxid, enableLowerCaseBases,
                    taxTree, dataNodes, fileNodes, idNodes, true, idStringGenerator, foldTaxaBelow);
            this.filter = filter;
        }

        @Override
        protected boolean handleStore(long kmer) {
            // Collect the contig's store value (matching what the DB fill would putLong: node.getTaxId()
            // when non-null). Guarded on a node change so it runs per contig, not per k-mer.
            if (node != lastCollectedNode) {
                lastCollectedNode = node;
                if (node != null && node.getTaxId() != null) {
                    values.add(node.getTaxId());
                }
            }
            // Lock-free combined membership-check-and-insert: the filter sets its bits atomically, so
            // the previous global lock on the filter is no longer needed. Each k-mer reported as new is
            // added to this reader's own radix-bucket counter; summed across readers these give the exact
            // number of distinct k-mers inserted.
            if (filter.putLong(kmer)) {
                bucketSizes[RadixKMerStore.radixOf(kmer, radixBits)]++;
                return true;
            }
            return false;
        }

        /**
         * Returns this reader's thread-local per-radix-bucket k-mer counts.
         *
         * @return the per-radix-bucket counts (length {@code 2^radixBits})
         */
        public int[] getBucketSizes() {
            return bucketSizes;
        }

        /**
         * Returns this reader's thread-local set of distinct store values (tax id strings).
         *
         * @return the distinct store values seen by this reader
         */
        public Set<String> getValues() {
            return values;
        }
    }
}