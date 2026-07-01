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
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.metagene.genestrip.ExecutionContext;
import org.metagene.genestrip.GSConfigKey;
import org.metagene.genestrip.GSGoalKey;
import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.bloom.BlockedKMerBloomFilter;
import org.metagene.genestrip.bloom.KMerProbFilter;
import org.metagene.genestrip.bloom.MurmurKMerBloomFilter;
import org.metagene.genestrip.bloom.XORKMerBloomFilter;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.refseq.AbstractRefSeqFastaReader;
import org.metagene.genestrip.refseq.AbstractStoreFastaReader;
import org.metagene.genestrip.refseq.AccessionMap;
import org.metagene.genestrip.refseq.RefSeqCategory;
import org.metagene.genestrip.store.RadixKMerStore;
import org.metagene.genestrip.tax.Rank;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;
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
         * Creates a DB size holder.
         *
         * @param size the total estimated number of distinct k-mers
         * @param bucketSizes the per-radix-bucket counts, or {@code null} if not computed
         */
        public DBSize(long size, int[] bucketSizes) {
            this.size = size;
            this.bucketSizes = bucketSizes;
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
    }

    private final ObjectGoal<AccessionMap, P> accessionMapGoal;
    private final ObjectGoal<Long, P> sizeGoal;

    private KMerProbFilter filter;
    // Number of low k-mer bits used as the radix (from config); the radix store created later must
    // use the same value. Set in doMakeThis().
    private int radixBits;
    // Deduplicated k-mer count per RadixKMerStore radix bucket (low radixBits bits of the k-mer),
    // computed alongside the total so a RadixKMerStore can later be sized per bucket. Shared by all
    // reader threads; see MyFastaReader.handleStore() for the thread-safe update.
    private int[] bucketSizes;

    private final boolean multiThreading;

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
     * @param sizeGoal the goal supplying the expected k-mer count for filter sizing
     * @param deps the additional goals this goal depends on
     */
    @SafeVarargs
    public FillBloomFilterGoal(P project, ExecutionContext bundle, ObjectGoal<Set<RefSeqCategory>, P> categoriesGoal,
                               ObjectGoal<Set<TaxIdNode>, P> taxNodesGoal, RefSeqFnaFilesDownloadGoal fnaFilesGoal,
                               ObjectGoal<Map<File, TaxIdNode>, P> additionalGoal,
                               ObjectGoal<AccessionMap, P> accessionMapGoal, ObjectGoal<Long, P> sizeGoal, Goal<P>... deps) {
        super(project, GSGoalKey.TEMPINDEX, bundle, categoriesGoal, taxNodesGoal, fnaFilesGoal, additionalGoal, Goal.append(deps, accessionMapGoal, sizeGoal));
        this.accessionMapGoal = accessionMapGoal;
        this.sizeGoal = sizeGoal;
        multiThreading = bundle.getThreads() > 0;
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
            if (tempFpp == BlockedKMerBloomFilter.DEFAULT_FPP) {
                filter = new BlockedKMerBloomFilter();
            }
            else {
                filter = booleanConfigValue(GSConfigKey.XOR_BLOOM_HASH) ?
                        new XORKMerBloomFilter(tempFpp) :
                        new MurmurKMerBloomFilter(tempFpp);
            }
            filter.ensureExpectedSize(sizeGoal.get(), false);
            radixBits = intConfigValue(GSConfigKey.RADIX_STORE_BITS);
            bucketSizes = new int[1 << radixBits];
            logHeapInfo();
            readFastas();
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
            for (int i = 0; i < bucketSizes.length; i++) {
                bucketSizes[i] = (int) (bucketSizes[i] / divisor);
            }
            set(new DBSize(entries, bucketSizes));
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
        return new MyFastaReader(intConfigValue(GSConfigKey.FASTA_LINE_SIZE_BYTES),
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
                booleanConfigValue(GSConfigKey.ENABLE_LOWERCASE_BASES));
    }

    /**
     * FASTA reader that adds each not-yet-seen k-mer to the temporary Bloom filter and counts the
     * distinct k-mers per radix bucket in a thread-safe way.
     */
    protected class MyFastaReader extends AbstractStoreFastaReader {
        private final KMerProbFilter filter;

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
                             KMerProbFilter filter, int maxGenomesPerTaxId, Rank maxGenomesPerTaxIdRank, long maxKmersPerTaxId, int maxDust, int stepSize, boolean completeGenomesOnly, StringLong2DigitTrie regionsPerTaxid, boolean enableLowerCaseBases) {
            super(bufferSize, taxNodes, accessionMap, k, maxGenomesPerTaxId, maxGenomesPerTaxIdRank, maxKmersPerTaxId, maxDust, stepSize, completeGenomesOnly, regionsPerTaxid, enableLowerCaseBases);
            this.filter = filter;
        }

        @Override
        protected boolean handleStore() {
            long kmer = byteRingBuffer.getStandardKMer();
            if (!filter.containsLong(kmer)) {
                int bucket = RadixKMerStore.radixOf(kmer, radixBits);
                if (multiThreading) {
                    synchronized (filter) {
                        // This is a trick to enable more parallelism -
                        // check again after synchronized to avoid synchronized further outside...
                        if (!filter.containsLong(kmer)) {
                            filter.putLong(kmer);
                            // Counted under the same lock as the filter update, so the per-bucket
                            // counts stay exactly in step with filter.getEntries() across threads
                            // (and race losers, whose re-check fails, are not counted).
                            bucketSizes[bucket]++;
                            return true;
                        }
                    }
                } else {
                    filter.putLong(kmer);
                    bucketSizes[bucket]++;
                    return true;
                }
            }
            return false;
        }
    }
}