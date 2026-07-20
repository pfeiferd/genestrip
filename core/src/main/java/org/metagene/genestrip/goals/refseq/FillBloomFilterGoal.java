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
import org.metagene.genestrip.bloom.CountingKMerProbFilter;
import org.metagene.genestrip.bloom.KMerProbFilter;
import org.metagene.genestrip.bloom.MurmurKMerBloomFilter;
import org.metagene.genestrip.bloom.XORKMerBloomFilter;
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

        /** The total estimated number of distinct k-mers. */
        private final long size;
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
         * @param size the total estimated number of distinct k-mers
         * @param bucketSizes the per-radix-bucket counts, or {@code null} if not computed
         * @param values the distinct fill values collected up front, or {@code null} if not collected
         */
        public DBSize(long size, int[] bucketSizes, Set<String> values) {
            this.size = size;
            this.bucketSizes = bucketSizes;
            this.values = values;
        }

        /**
         * Returns the total estimated number of distinct k-mers.
         *
         * @return the total estimated number of distinct k-mers
         */
        public long getSize() {
            return size;
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
    private final ObjectGoal<Long, P> sizeGoal;

    // The temporary size-estimation filter is always an XOR/Murmur filter (both count their
    // entries), and its entry count is read below to derive the store sizing.
    private CountingKMerProbFilter filter;
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
     * @param sizeGoal the goal supplying the expected k-mer count for filter sizing
     * @param deps the additional goals this goal depends on
     */
    @SafeVarargs
    public FillBloomFilterGoal(P project, ExecutionContext bundle, ObjectGoal<Set<RefSeqCategory>, P> categoriesGoal,
                               ObjectGoal<Set<TaxIdNode>, P> taxNodesGoal, RefSeqFnaFilesDownloadGoal fnaFilesGoal,
                               ObjectGoal<Map<File, TaxIdNode>, P> additionalGoal,
                               ObjectGoal<AccessionMap, P> accessionMapGoal, ObjectGoal<TaxTree, P> taxTreeGoal,
                               ObjectGoal<Long, P> sizeGoal, Goal<P>... deps) {
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
        try {
            double tempFpp = doubleConfigValue(GSConfigKey.TEMP_BLOOM_FILTER_FPP);
            // The temporary size-estimation filter is filled concurrently by the reader threads, so it
            // must support a thread-safe putLongIfAbsent; XOR/Murmur do (via their bit vector's bucket
            // locks). BlockedKMerBloomFilter is deliberately not used here.
            filter = booleanConfigValue(GSConfigKey.XOR_BLOOM_HASH) ?
                    new XORKMerBloomFilter(tempFpp, sizeGoal.get()) :
                    new MurmurKMerBloomFilter(tempFpp, sizeGoal.get());
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
            long entries = filter.getEntries();
            // We have to account for the missing entries in the bloom filter due to
            // inherent FPP. The formula from below works really well,
            // so we can allow for a low FPP for the bloom filter from 'bloomFilterGoal'
            // and save memory during db construction.
            // It is a conservative estimate too, since collisions occur in the process
            // of filling (as opposed to the FPP formula that considers a filled
            // bloom filter).
            double divisor = 1d - doubleConfigValue(GSConfigKey.TEMP_BLOOM_FILTER_FPP);
            entries = (long) (entries / divisor);
            // Apply the same FPP correction per bucket so the bucket sizes stay consistent with the
            // total (their sum stays approximately 'entries').
            int[] correctedBucketSizes = new int[bucketSizes.length];
            for (int i = 0; i < correctedBucketSizes.length; i++) {
                correctedBucketSizes[i] = (int) (bucketSizes[i] / divisor);
            }
            set(new DBSize(entries, correctedBucketSizes, collectedValues));
        } catch (IOException e) {
            throw new RuntimeException(e);
        } finally {
            filter = null;
            System.gc(); // Time to run GC after freeing up the bloom filter's memory.
            cleanUpThreads();
        }
    }

    @Override
    protected void afterReadFastas(AbstractRefSeqFastaReader.StringLong2DigitTrie regionsPerTaxid) {
        if (getLogger().isInfoEnabled()) {
            getLogger().info("Bloom filter size in kmers: " + filter.getEntries());
            getLogger().info("Duplication factor: " + ((double) sizeGoal.get()) / filter.getEntries());
        }
        if (getLogger().isDebugEnabled()) {
            List<StringLongDigitTrie.StringLong> list = new ArrayList<>();
            regionsPerTaxid.collect(list);
            getLogger().debug("Regions ber taxid:");
            getLogger().debug(list);
        }
    }

    @Override
    protected AbstractStoreFastaReader createFastaReader(AbstractRefSeqFastaReader.StringLong2DigitTrie regionsPerTaxid) {
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
                taxNodesGoal.get(),
                isIncludeRefSeqFna() ? accessionMapGoal.get() : null,
                intConfigValue(GSConfigKey.KMER_SIZE),
                filter,
                intConfigValue(GSConfigKey.MAX_GENOMES_PER_TAXID),
                (Rank) configValue(GSConfigKey.MAX_GENOMES_PER_TAXID_RANK),
                longConfigValue(GSConfigKey.MAX_KMERS_PER_TAXID),
                intConfigValue(GSConfigKey.MAX_DUST),
                intConfigValue(GSConfigKey.STEP_SIZE),
                booleanConfigValue(GSConfigKey.COMPLETE_GENOMES_ONLY),
                regionsPerTaxid,
                booleanConfigValue(GSConfigKey.ENABLE_LOWERCASE_BASES),
                taxTreeGoal.get(),
                booleanConfigValue(GSConfigKey.DATA_NODES),
                booleanConfigValue(GSConfigKey.FILE_NODES),
                booleanConfigValue(GSConfigKey.ID_NODES),
                idStringGenerator);
        readers.add(fastaReader);
        return fastaReader;
    }

    /**
     * FASTA reader that adds each not-yet-seen k-mer to the temporary Bloom filter and counts the
     * distinct k-mers per radix bucket in a thread-safe way.
     */
    protected class MyFastaReader extends ReworkingStoreFastaReader {
        private final KMerProbFilter filter;
        // Thread-local per-radix-bucket k-mer counts and distinct store values - this reader runs on a
        // single thread, so no synchronization is needed; doMakeThis() merges all readers afterwards.
        private final int[] bucketSizes = new int[1 << radixBits];
        private final Set<String> values = new HashSet<>();
        // Last node whose taxid was collected, so the (rare) region-boundary collection skips the
        // repeated per-k-mer adds within a region.
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
         * @param maxGenomesPerTaxId the maximum number of genomes kept per tax id
         * @param maxGenomesPerTaxIdRank the rank at which the genome limit is applied
         * @param maxKmersPerTaxId the maximum number of k-mers kept per tax id
         * @param maxDust the maximum allowed low-complexity (dust) run length
         * @param stepSize the k-mer sampling step size
         * @param completeGenomesOnly whether only complete genomes are considered
         * @param regionsPerTaxid the per-tax-id region counter
         * @param enableLowerCaseBases whether lower-case bases are processed
         */
        public MyFastaReader(int bufferSize, Set<TaxIdNode> taxNodes, AccessionMap accessionMap, int k,
                             KMerProbFilter filter, int maxGenomesPerTaxId, Rank maxGenomesPerTaxIdRank, long maxKmersPerTaxId, int maxDust, int stepSize, boolean completeGenomesOnly, StringLong2DigitTrie regionsPerTaxid, boolean enableLowerCaseBases,
                             TaxTree taxTree, boolean dataNodes, boolean fileNodes, boolean idNodes, IDStringGenerator idStringGenerator) {
            super(bufferSize, taxNodes, accessionMap, k, maxGenomesPerTaxId, maxGenomesPerTaxIdRank, maxKmersPerTaxId, maxDust, stepSize, completeGenomesOnly, regionsPerTaxid, enableLowerCaseBases,
                    taxTree, dataNodes, fileNodes, idNodes, true, idStringGenerator);
            this.filter = filter;
        }

        @Override
        protected boolean handleStore() {
            // Collect the region's store value (matching what the DB fill would putLong: node.getTaxId()
            // when non-null). Guarded on a node change so it runs per region, not per k-mer.
            if (node != lastCollectedNode) {
                lastCollectedNode = node;
                if (node != null && node.getTaxId() != null) {
                    values.add(node.getTaxId());
                }
            }
            long kmer = byteRingBuffer.getStandardKMer();
            // Lock-free combined membership-check-and-insert: the filter sets its bits atomically, so
            // the previous global lock on the filter is no longer needed. Each k-mer reported as new is
            // added to this reader's own radix-bucket counter; summed across readers these give the exact
            // number of distinct k-mers inserted (whereas filter.getEntries() counts every insert call,
            // i.e. k-mer occurrences including duplicates).
            if (filter.putLongIfAbsent(kmer)) {
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