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
package org.metagene.genestrip.finertree.goals;

import it.unimi.dsi.fastutil.objects.Object2LongMap;
import it.unimi.dsi.fastutil.objects.Object2LongOpenHashMap;
import org.metagene.genestrip.ExecutionContext;
import org.metagene.genestrip.GSConfigKey;
import org.metagene.genestrip.finertree.FTConfigKey;
import org.metagene.genestrip.finertree.FTGoalKey;
import org.metagene.genestrip.finertree.FTProject;
import org.metagene.genestrip.finertree.bloom.XORKMerIndexBloomFilter;
import org.metagene.genestrip.goals.refseq.FastaReaderGoal;
import org.metagene.genestrip.goals.refseq.RefSeqFnaFilesDownloadGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.refseq.AbstractRefSeqFastaReader;
import org.metagene.genestrip.refseq.AbstractStoreFastaReader;
import org.metagene.genestrip.refseq.AccessionMap;
import org.metagene.genestrip.refseq.RefSeqCategory;
import org.metagene.genestrip.store.Database;
import org.metagene.genestrip.store.KMerSortedArray;
import org.metagene.genestrip.tax.Rank;
import org.metagene.genestrip.tax.SmallTaxTree;
import org.metagene.genestrip.tax.TaxTree;

import java.io.File;
import java.io.IOException;
import java.util.*;

public class KMersPerDBTaxidGoal<P extends FTProject> extends FastaReaderGoal<Object2LongMap<String>, P> implements Goal.LogHeapInfo {
    private final ObjectGoal<AccessionMap, P> accessionMapGoal;
    private final ObjectGoal<Database, P> storeGoal;

    private KMerSortedArray<String> kMerSortedArray;
    private long[] counts;
    private XORKMerIndexBloomFilter filter;

    @SafeVarargs
    public KMersPerDBTaxidGoal(P project, FTGoalKey key, ExecutionContext bundle, ObjectGoal<Set<RefSeqCategory>, P> categoriesGoal,
                              ObjectGoal<Set<TaxTree.TaxIdNode>, P> taxNodesGoal,
                              RefSeqFnaFilesDownloadGoal fnaFilesGoal,
                              ObjectGoal<Map<File, TaxTree.TaxIdNode>, P> additionalGoal,
                              ObjectGoal<AccessionMap, P> accessionMapGoal, ObjectGoal<Database, P> storeGoal,
                              Goal<P>... deps) {
        super(project, key, bundle, categoriesGoal, taxNodesGoal, fnaFilesGoal, additionalGoal, Goal.append(deps, accessionMapGoal, storeGoal));
        this.storeGoal = storeGoal;
        this.accessionMapGoal = accessionMapGoal;
    }

    @Override
    protected void doMakeThis() {
        SmallTaxTree tree = storeGoal.get().getTaxTree();
        Object2LongMap stats = storeGoal.get().getStats();
        // Estimate the filter size by summing up from leaves to root for each leaf taxid in the DB.
        long sum = 0;
        Iterator<SmallTaxTree.SmallTaxIdNode> iterator = tree.iterator();
        while (iterator.hasNext()) {
            SmallTaxTree.SmallTaxIdNode node = iterator.next();
            // Leaves where originally stored...
            if (FTConsistencyGoal.acceptNode(node)) {
                sum += FTConsistencyGoal.getPathSum(node, stats);
            }
        }
        filter = new XORKMerIndexBloomFilter(doubleConfigValue(FTConfigKey.FT_BLOOM_FILTER_FPP));
        long bitSize = filter.ensureExpectedSize(sum, false);
        if (getLogger().isInfoEnabled()) {
            getLogger().info("Filter size in MB: " + (bitSize / 8 / 1024 / 1024));
        }
        kMerSortedArray = storeGoal.get().getKmerStore();
        counts = new long[kMerSortedArray.getNValues()];

        try {
            readFastas();
        } catch (IOException e) {
            throw new RuntimeException(e);
        } finally {
            cleanUpThreads();
        }

        Object2LongOpenHashMap<String> map = new Object2LongOpenHashMap<>(counts.length);
        for (int i = 0; i < counts.length; i++) {
            map.put(kMerSortedArray.getValueForIndex(i), counts[i]);
        }
        set(map);
    }

    @Override
    protected AbstractStoreFastaReader createFastaReader(AbstractRefSeqFastaReader.StringLong2DigitTrie regionsPerTaxid) {
        return new KMersPerDBTaxidGoal.MyFastaReader(intConfigValue(GSConfigKey.FASTA_LINE_SIZE_BYTES),
                taxNodesGoal.get(),
                isIncludeRefSeqFna() ? accessionMapGoal.get() : null,
                intConfigValue(GSConfigKey.KMER_SIZE),
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
        public MyFastaReader(int bufferSize, Set<TaxTree.TaxIdNode> taxNodes, AccessionMap accessionMap,
                             int k, int maxGenomesPerTaxId, Rank maxGenomesPerTaxIdRank, long maxKmersPerTaxId, int maxDust, int stepSize, boolean completeGenomesOnly, StringLong2DigitTrie regionsPerTaxid, boolean enableLowerCaseBases) {
            super(bufferSize, taxNodes, accessionMap, k, maxGenomesPerTaxId, maxGenomesPerTaxIdRank, maxKmersPerTaxId, maxDust, stepSize, completeGenomesOnly, regionsPerTaxid, enableLowerCaseBases);
        }

        @Override
        protected boolean handleStore() {
            long kmer = byteRingBuffer.getStandardKMer();
            String taxid = kMerSortedArray.getLong(kmer, null);
            // Checks whether k-mer is stored at all
            if (taxid != null) {
                // We only consider kmer above the node's taxid.
                if (!node.getTaxId().equals(taxid)) {
                    int index = kMerSortedArray.getIndexForValue(node.getTaxId());
                    // Checks whether it's a duplicate under that taxid.
                    if (!filter.containsLongInt(kmer, index)) {
                        synchronized (counts) {
                            filter.putLongInt(kmer, index);
                            counts[index]++;
                        }
                    }
                    return true;
                }
            }
            return false;
        }
    }
}
