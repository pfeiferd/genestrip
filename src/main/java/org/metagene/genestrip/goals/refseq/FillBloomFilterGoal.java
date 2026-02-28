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
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.metagene.genestrip.ExecutionContext;
import org.metagene.genestrip.GSConfigKey;
import org.metagene.genestrip.GSGoalKey;
import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.bloom.KMerProbFilter;
import org.metagene.genestrip.bloom.MurmurKMerBloomFilter;
import org.metagene.genestrip.bloom.XORKMerBloomFilter;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.refseq.AbstractRefSeqFastaReader;
import org.metagene.genestrip.refseq.AbstractStoreFastaReader;
import org.metagene.genestrip.refseq.AccessionMap;
import org.metagene.genestrip.refseq.RefSeqCategory;
import org.metagene.genestrip.tax.Rank;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;
import org.metagene.genestrip.util.StringLongDigitTrie;

// Not needed anymore due to HyperLogLog estimation of DB size.
//@Deprecated
public class FillBloomFilterGoal<P extends GSProject> extends FastaReaderGoal<Long, P> implements Goal.LogHeapInfo {
    private final ObjectGoal<AccessionMap, P> accessionMapGoal;
    private final ObjectGoal<Long, P> sizeGoal;

    private KMerProbFilter filter;

    private final boolean multiThreading;

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

    @Override
    protected void allDependentsMade() {
        // To save memory...
        doCleanThis();
    }

    @Override
    protected void doMakeThis() {
        try {
            filter = booleanConfigValue(GSConfigKey.XOR_BLOOM_HASH) ?
                    new XORKMerBloomFilter(doubleConfigValue(GSConfigKey.TEMP_BLOOM_FILTER_FPP)) :
                    new MurmurKMerBloomFilter(doubleConfigValue(GSConfigKey.TEMP_BLOOM_FILTER_FPP));
            filter.ensureExpectedSize(sizeGoal.get(), false);
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
            entries = (long) (entries / (1d - doubleConfigValue(GSConfigKey.TEMP_BLOOM_FILTER_FPP)));
            set(entries);
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

    protected class MyFastaReader extends AbstractStoreFastaReader {
        private final KMerProbFilter filter;

        public MyFastaReader(int bufferSize, Set<TaxIdNode> taxNodes, AccessionMap accessionMap, int k,
                             KMerProbFilter filter, int maxGenomesPerTaxId, Rank maxGenomesPerTaxIdRank, long maxKmersPerTaxId, int maxDust, int stepSize, boolean completeGenomesOnly, StringLong2DigitTrie regionsPerTaxid, boolean enableLowerCaseBases) {
            super(bufferSize, taxNodes, accessionMap, k, maxGenomesPerTaxId, maxGenomesPerTaxIdRank, maxKmersPerTaxId, maxDust, stepSize, completeGenomesOnly, regionsPerTaxid, enableLowerCaseBases);
            this.filter = filter;
        }

        @Override
        protected boolean handleStore() {
            long kmer = byteRingBuffer.getStandardKMer();
            if (!filter.containsLong(kmer)) {
                if (multiThreading) {
                    synchronized (filter) {
                        // This is a trick to enable more parallelism -
                        // check again after synchronized to avoid synchronized further outside...
                        if (!filter.containsLong(kmer)) {
                            filter.putLong(kmer);
                            return true;
                        }
                    }
                } else {
                    filter.putLong(kmer);
                    return true;
                }
            }
            return false;
        }
    }
}