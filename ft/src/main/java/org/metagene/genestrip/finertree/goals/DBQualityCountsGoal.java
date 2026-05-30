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
import org.metagene.genestrip.ExecutionContext;
import org.metagene.genestrip.GSConfigKey;
import org.metagene.genestrip.finertree.FTConfigKey;
import org.metagene.genestrip.finertree.FTGoalKey;
import org.metagene.genestrip.finertree.FTProject;
import org.metagene.genestrip.finertree.bloom.XORKMerIndexBloomFilter;
import org.metagene.genestrip.finertree.refseq.AbstractUpdateFastaReader;
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
import java.io.Serializable;
import java.util.*;

public class DBQualityCountsGoal<P extends FTProject> extends FastaReaderGoal<Map<String, DBQualityCountsGoal.Counts>, P> implements Goal.LogHeapInfo {
    private final ObjectGoal<AccessionMap, P> accessionMapGoal;
    private final ObjectGoal<Database, P> storeGoal;

    private SmallTaxTree tree;
    private KMerSortedArray<SmallTaxTree.SmallTaxIdNode> kMerSortedArray;
    private XORKMerIndexBloomFilter filter;
    private Map<String, Counts> map;

    @SafeVarargs
    public DBQualityCountsGoal(P project, FTGoalKey key, ExecutionContext bundle, ObjectGoal<Set<RefSeqCategory>, P> categoriesGoal,
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
        try {
            map = new HashMap<>();
            tree = storeGoal.get().getTaxTree();
            Object2LongMap<String> stats = storeGoal.get().getStats();
            // Estimate the filter size by summing up from species to root for each species in the DB.
            // It is a highly conservative estimate because k-mers on ranks above species are hardly
            // ever shared more than thrice (as found by measuring).
            long size = 0;
            Iterator<SmallTaxTree.SmallTaxIdNode> iterator = tree.iterator();
            while (iterator.hasNext()) {
                SmallTaxTree.SmallTaxIdNode node = iterator.next();
                if (node.getSubNodes() == null) {
                    size += getPathSum(node, stats);
                }
            }
            filter = new XORKMerIndexBloomFilter(doubleConfigValue(FTConfigKey.FT_BLOOM_FILTER_FPP));
            long bitSize = filter.ensureExpectedSize(size, false);
            if (getLogger().isInfoEnabled()) {
                getLogger().info("Filter size in MB: " + (bitSize / 8 / 1024 / 1024));
            }
            kMerSortedArray = storeGoal.get().convertKMerStore();

            readFastas();
            if (getLogger().isInfoEnabled()) {
                getLogger().info("Actual filter entries: " + filter.getEntries());
            }
            if (getLogger().isErrorEnabled()) {
                if (filter.getEntries() > 2 * size) {
                    getLogger().error("Entries exceed filter size by over factor 2. Something went wrong!");
                }
            }
            // Count tp plus fp
            // Add k-mers from species upwards for each species:
            for (String taxid : map.keySet()) {
                Counts counts = map.get(taxid);
                SmallTaxTree.SmallTaxIdNode node = tree.getNodeByTaxId(taxid);
                if (node != null) {
                    counts.tpPlusFp += getPathSum(node, stats);
                } else if (getLogger().isWarnEnabled()) {
                    getLogger().warn("Missing node for taxid " + taxid + ".");
                }
            }

            // We need a separate map for averaging so that aggregations do not get mixed up e.g. between genus and species.
            Map<String, Counts> aggMap = new HashMap<>();
            List<Rank> aggRanks = Arrays.asList(Rank.CELLULAR_ROOT, Rank.ACELLULAR_ROOT, Rank.SPECIES, Rank.GENUS);
            iterator = tree.iterator();
            while (iterator.hasNext()) {
                SmallTaxTree.SmallTaxIdNode node = iterator.next();
                Counts counts = map.get(node.getTaxId());
                if (counts != null) {
                    for (Rank rank : aggRanks) {
                        SmallTaxTree.SmallTaxIdNode rankedNode = toRankedNode(node, rank);
                        // No aggregation for node who are already contained in map.
                        if (rankedNode != null && !map.containsKey(rankedNode.getTaxId())) {
                            Counts c = aggMap.get(rankedNode.getTaxId());
                            if (c == null) {
                                c = new Counts();
                                aggMap.put(rankedNode.getTaxId(), c);
                            }
                            // Leads to a weighted average:
                            // Node weight is proportional to positives of each aggregated node
                            c.aggregate(counts);
                        }
                    }
                }
            }
            map.putAll(aggMap);
            set(map);
        } catch (IOException e) {
            throw new RuntimeException(e);
        } finally {
            cleanUpThreads();
        }
    }

    @Override
    protected AbstractStoreFastaReader createFastaReader(AbstractRefSeqFastaReader.StringLong2DigitTrie regionsPerTaxid) {
        return new MyFastaReader(intConfigValue(GSConfigKey.FASTA_LINE_SIZE_BYTES),
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

    protected class MyFastaReader extends AbstractUpdateFastaReader {
        public MyFastaReader(int bufferSize, Set<TaxTree.TaxIdNode> taxNodes, AccessionMap accessionMap,
                             int k, int maxGenomesPerTaxId, Rank maxGenomesPerTaxIdRank, long maxKmersPerTaxId, int maxDust, int stepSize, boolean completeGenomesOnly, StringLong2DigitTrie regionsPerTaxid, boolean enableLowerCaseBases) {
            super(bufferSize, taxNodes, accessionMap, k, maxGenomesPerTaxId, maxGenomesPerTaxIdRank, maxKmersPerTaxId, maxDust, stepSize, completeGenomesOnly, regionsPerTaxid, enableLowerCaseBases, booleanConfigValue(GSConfigKey.ID_NODES), booleanConfigValue(GSConfigKey.FILE_NODES), booleanConfigValue(GSConfigKey.DATA_NODES));
        }

        @Override
        protected SmallTaxTree getTree() {
            return tree;
        }

        @Override
        protected boolean handleStore() {
            // There may be no corresponding node in the database:
            if (leafNode != null) {
                long kmer = byteRingBuffer.getStandardKMer();
                SmallTaxTree.SmallTaxIdNode storedNode = kMerSortedArray.getLong(kmer, null);
                // Checks whether k-mer is stored in database at all:
                if (storedNode != null) {
                    int index = kMerSortedArray.getIndexForValue(leafNode);
                    if (index >= 0) {
                        synchronized (filter) {
                            // Checks whether it's a duplicate under that taxid.
                            if (!filter.containsLongInt(kmer, index)) {
                                filter.putLongInt(kmer, index);
                                Counts counts = map.get(leafNode.getTaxId());
                                if (counts == null) {
                                    counts = new Counts();
                                    map.put(leafNode.getTaxId(), counts);
                                }
                                counts.tpPlusFn++;
                                // Is the stored node on path of the file's node?
                                // Moving up from the file node to the potentially higher ranked stored node.
                                SmallTaxTree.SmallTaxIdNode pathNode = leafNode;
                                while (pathNode != null && pathNode != storedNode) {
                                    pathNode = pathNode.getParent();
                                }
                                if (pathNode == storedNode) {
                                    // Stored node on path from file node: True positive
                                    counts.tp++;
                                }
                                return true;
                            }
                        }
                    } else if (getLogger().isWarnEnabled()) {
                        getLogger().warn("No kmer-index for taxid " + storedNode.getTaxId() + " found.");
                    }
                }
            }
            return false;
        }
    }

    protected SmallTaxTree.SmallTaxIdNode toRankedNode(SmallTaxTree.SmallTaxIdNode node, Rank r) {
        while (node != null && !r.equals(node.getRank())) {
            node = node.getParent();
        }
        return node;
    }

    protected long getPathSum(SmallTaxTree.SmallTaxIdNode node, Object2LongMap stats) {
        long res = 0L;
        for (; node != null; node = node.getParent()) {
            res += stats.getOrDefault(node.getTaxId(), 0L);
        }
        return res;
    }

    public static class Counts implements Serializable {
        private static final long serialVersionUID = 1L;

        private long tpPlusFp;
        private long tp;
        private long tpPlusFn;
        private int aggregations;
        private double aggPrecisionSum;
        private double aggRecallSum;

        public Counts() {
        }

        public long getTpPlusFn() {
            return tpPlusFn;
        }

        public long getTp() {
            return tp;
        }

        public long getTpPlusFp() {
            return tpPlusFp;
        }

        // Weighted
        public double getPrecision() {
            return ((double) tp / tpPlusFp);
        }

        // Unweighted
        public double getAvgPrecision() {
            if (aggPrecisionSum == 0) {
                return getPrecision();
            }
            else {
                return aggPrecisionSum / aggregations;
            }
        }

        // Weighted
        public double getRecall() {
            return ((double) tp / tpPlusFn);
        }

        // Unweighted
        public double getAvgRecall() {
            if (aggPrecisionSum == 0) {
                return getRecall();
            }
            else {
                return aggRecallSum / aggregations;
            }
        }

        private void aggregate(Counts counts) {
            tp += counts.tp;
            tpPlusFp += counts.tpPlusFp;
            tpPlusFn += counts.tpPlusFn;
            aggregations++;
            aggPrecisionSum += counts.getAvgPrecision();
            aggRecallSum += counts.getAvgRecall();
        }
    }
}
