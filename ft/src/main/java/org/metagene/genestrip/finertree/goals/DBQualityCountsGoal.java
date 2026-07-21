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
import org.metagene.genestrip.store.KMerStore;
import org.metagene.genestrip.tax.Rank;
import org.metagene.genestrip.tax.SmallTaxTree;
import org.metagene.genestrip.tax.TaxTree;

import java.io.File;
import java.io.IOException;
import java.io.Serializable;
import java.util.*;

/**
 * Computes intrinsic database-quality counts per tax id by re-reading the underlying genomic fasta
 * files and comparing the *k*-mers they contain against those stored in the database. For each tax id
 * it accumulates true positives, true-positives-plus-false-positives and true-positives-plus-false-
 * negatives (from which precision and recall are derived), aggregating results up selected ranks. A
 * XOR bloom filter is used to detect duplicate (k-mer, tax id) pairs.
 *
 * @param <P> the concrete FT project type
 */
public class DBQualityCountsGoal<P extends FTProject> extends FastaReaderGoal<Map<String, DBQualityCountsGoal.Counts>, P> implements Goal.LogHeapInfo {
    private final ObjectGoal<AccessionMap, P> accessionMapGoal;
    private final ObjectGoal<Database, P> storeGoal;

    private SmallTaxTree tree;
    private KMerStore<SmallTaxTree.SmallTaxIdNode> kMerSortedArray;
    private XORKMerIndexBloomFilter filter;
    private Map<String, Counts> map;
    private long entries;

    /**
     * Creates the goal, depending on the accession map and the loaded database in addition to the
     * standard fasta-reader dependencies.
     *
     * @param project          the FT project
     * @param key              the goal key
     * @param bundle           the execution context used to run the fasta readers
     * @param categoriesGoal   the goal providing the RefSeq categories to read
     * @param taxNodesGoal     the goal providing the tax nodes to be included
     * @param fnaFilesGoal     the goal providing the downloaded genomic fasta files
     * @param additionalGoal   the goal providing additional fasta files mapped to tax nodes
     * @param accessionMapGoal the goal providing the accession-to-tax-node map
     * @param storeGoal        the goal providing the loaded database
     * @param deps             further goals this goal depends on
     */
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

    /**
     * Re-reads the genomic fasta files, comparing their *k*-mers against the database to accumulate the
     * per-tax-id true-positive and positive counts, aggregates the counts up the selected ranks and
     * stores the resulting map as this goal's value.
     */
    @Override
    protected void doMakeThis() {
        try {
            map = new HashMap<>();
            tree = storeGoal.get().getTaxTree();
            entries = 0;
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
            filter = new XORKMerIndexBloomFilter(doubleConfigValue(FTConfigKey.FT_BLOOM_FILTER_FPP), size);
            long bitSize = filter.getBitSize();
            if (getLogger().isInfoEnabled()) {
                getLogger().info("Filter size in MB: " + (bitSize / 8 / 1024 / 1024));
            }
            kMerSortedArray = storeGoal.get().convertKMerStore();

            readFastas();
            if (getLogger().isInfoEnabled()) {
                getLogger().info("Filter entries: " + entries);
            }
            if (getLogger().isErrorEnabled()) {
                if (entries > 2 * size) {
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
            map = null;
            tree = null;
            kMerSortedArray = null;
            filter = null;
            cleanUpThreads();
        }
    }

    /**
     * Creates the fasta reader that compares genome *k*-mers against the database, configured from the
     * project's configuration values.
     *
     * @param regionsPerTaxid the trie counting regions per tax id
     * @return the fasta reader to use for reading the genomic fasta files
     */
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

    /**
     * Fasta reader that, for each *k*-mer read from a genome, checks whether it is stored in the
     * database and whether the stored node lies on the path from the read's leaf node, updating the
     * per-tax-id counts accordingly while deduplicating via the bloom filter.
     */
    protected class MyFastaReader extends AbstractUpdateFastaReader {
        /**
         * Creates the reader, taking the id/file/data node flags from the goal's configuration.
         *
         * @param bufferSize             the input read buffer size
         * @param taxNodes               the tax nodes to be included
         * @param accessionMap           the map from accession numbers to tax nodes
         * @param k                      the k-mer length
         * @param maxGenomesPerTaxId     the maximum number of genomes to consider per tax id
         * @param maxGenomesPerTaxIdRank the rank at which the per-tax-id genome limit applies
         * @param maxKmersPerTaxId       the maximum number of k-mers to store per tax id
         * @param maxDust                the maximum dust (low-complexity) threshold
         * @param stepSize               the step size between stored k-mers
         * @param completeGenomesOnly    whether only complete genomes are considered
         * @param regionsPerTaxid        the trie counting regions per tax id
         * @param enableLowerCaseBases   whether lower-case bases are treated as regular bases
         */
        public MyFastaReader(int bufferSize, Set<TaxTree.TaxIdNode> taxNodes, AccessionMap accessionMap,
                             int k, int maxGenomesPerTaxId, Rank maxGenomesPerTaxIdRank, long maxKmersPerTaxId, int maxDust, int stepSize, boolean completeGenomesOnly, StringLong2DigitTrie regionsPerTaxid, boolean enableLowerCaseBases) {
            super(bufferSize, taxNodes, accessionMap, k, maxGenomesPerTaxId, maxGenomesPerTaxIdRank, maxKmersPerTaxId, maxDust, stepSize, completeGenomesOnly, regionsPerTaxid, enableLowerCaseBases, booleanConfigValue(GSConfigKey.ID_NODES), booleanConfigValue(GSConfigKey.FILE_NODES), booleanConfigValue(GSConfigKey.DATA_NODES));
        }

        /**
         * Returns the taxonomy tree of the loaded database.
         *
         * @return the taxonomy tree
         */
        @Override
        protected SmallTaxTree getTree() {
            return tree;
        }

        /**
         * Looks up the current *k*-mer in the database and, if it is stored and has not already been
         * seen for the read's leaf node, records it in the bloom filter and updates the per-tax-id
         * counts (incrementing the true-positive count when the stored node lies on the path from the
         * leaf node).
         *
         * @return {@code true} if the *k*-mer was counted, {@code false} otherwise
         */
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
                        // Checks whether it's a duplicate under that taxid.
                        if (filter.putLongInt(kmer, index)) {
                            synchronized (map) {
                                entries++;
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
                            }
                            return true;
                        }
                    } else if (getLogger().isWarnEnabled()) {
                        getLogger().warn("No kmer-index for taxid " + storedNode.getTaxId() + " found.");
                    }
                }
            }
            return false;
        }
    }

    /**
     * Walks up the taxonomy tree until a node of the given rank is reached.
     *
     * @param node the node to start from
     * @param r    the rank to look for
     * @return the nearest ancestor of {@code node} (or the node itself) having rank {@code r}, or
     * {@code null} if there is none
     */
    protected SmallTaxTree.SmallTaxIdNode toRankedNode(SmallTaxTree.SmallTaxIdNode node, Rank r) {
        while (node != null && !r.equals(node.getRank())) {
            node = node.getParent();
        }
        return node;
    }

    /**
     * Sums the per-tax-id stored k-mer counts along the path from a node up to the root.
     *
     * @param node  the node to start from
     * @param stats the per-tax-id stored k-mer counts keyed by tax id
     * @return the sum of the per-tax-id stored *k*-mer counts from {@code stats} along the path from
     * {@code node} up to the root
     */
    protected long getPathSum(SmallTaxTree.SmallTaxIdNode node, Object2LongMap stats) {
        long res = 0L;
        for (; node != null; node = node.getParent()) {
            res += stats.getOrDefault(node.getTaxId(), 0L);
        }
        return res;
    }

    /**
     * Per-tax-id tally of true positives, true-positives-plus-false-positives and true-positives-plus-
     * false-negatives, plus aggregated precision/recall sums used to compute weighted and unweighted
     * averages.
     */
    public static class Counts implements Serializable {
        private static final long serialVersionUID = 1L;

        /**
         * True positives plus false positives (all k-mers stored under this tax id).
         */
        private long tpPlusFp;
        /**
         * True positives (k-mers correctly stored under this tax id).
         */
        private long tp;
        /**
         * True positives plus false negatives (all k-mers read from this tax id's genomes).
         */
        private long tpPlusFn;
        /**
         * Number of child nodes aggregated into this tally.
         */
        private int aggregations;
        /**
         * Sum of the aggregated children's precision values.
         */
        private double aggPrecisionSum;
        /**
         * Sum of the aggregated children's recall values.
         */
        private double aggRecallSum;

        /**
         * Creates an empty tally with all counts set to zero.
         */
        public Counts() {
        }

        /**
         * Returns the true-positives-plus-false-negatives count.
         *
         * @return the true positives plus false negatives
         */
        public long getTpPlusFn() {
            return tpPlusFn;
        }

        /**
         * Returns the true-positives count.
         *
         * @return the true positives
         */
        public long getTp() {
            return tp;
        }

        /**
         * Returns the true-positives-plus-false-positives count.
         *
         * @return the true positives plus false positives
         */
        public long getTpPlusFp() {
            return tpPlusFp;
        }

        // Weighted

        /**
         * Returns the weighted precision.
         *
         * @return the precision {@code tp / (tp + fp)}
         */
        public double getPrecision() {
            return ((double) tp / tpPlusFp);
        }

        // Unweighted

        /**
         * Returns the unweighted average precision.
         *
         * @return the unweighted average precision over aggregated nodes, or {@link #getPrecision()}
         * if there was no aggregation
         */
        public double getAvgPrecision() {
            if (aggPrecisionSum == 0) {
                return getPrecision();
            } else {
                return aggPrecisionSum / aggregations;
            }
        }

        // Weighted

        /**
         * Returns the weighted recall.
         *
         * @return the recall {@code tp / (tp + fn)}
         */
        public double getRecall() {
            return ((double) tp / tpPlusFn);
        }

        // Unweighted

        /**
         * Returns the unweighted average recall.
         *
         * @return the unweighted average recall over aggregated nodes, or {@link #getRecall()} if
         * there was no aggregation
         */
        public double getAvgRecall() {
            if (aggRecallSum == 0) {
                return getRecall();
            } else {
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
