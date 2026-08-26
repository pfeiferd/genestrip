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
import org.metagene.genestrip.finertree.probfilter.KMerIndexFilterHelper;
import org.metagene.genestrip.goals.refseq.RefSeqFnaFilesDownloadGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.probfilter.BlockedBloomFilter;
import org.metagene.genestrip.probfilter.ProbFilter;
import org.metagene.genestrip.refseq.AbstractRefSeqFastaReader;
import org.metagene.genestrip.refseq.AccessionMap;
import org.metagene.genestrip.refseq.RefSeqCategory;
import org.metagene.genestrip.store.Database;
import org.metagene.genestrip.tax.SmallTaxTree;
import org.metagene.genestrip.tax.TaxTree;
import org.metagene.genestrip.tax.TaxNodeSelection;

import java.io.File;
import java.io.IOException;
import java.io.Serializable;
import java.util.*;

/**
 * Computes intrinsic database-quality counts per tax id by re-reading the underlying genomic fasta
 * files and comparing the *k*-mers they contain against those stored in the database. For each tax id
 * it accumulates true positives, true-positives-plus-false-positives and true-positives-plus-false-
 * negatives (from which precision and recall are derived), aggregating results up selected ranks. A
 * bloom filter is used to detect duplicate (k-mer, leaf) pairs.
 * <p>
 * The reading itself lives in {@link AbstractDBQualityGoal}; the tallies are here, because only this
 * goal has any.
 *
 * @param <P> the concrete FT project type
 */
public class DBQualityCountsGoal<P extends FTProject> extends AbstractDBQualityGoal<Map<String, DBQualityCountsGoal.Counts>, P> {
    /** Supplies the number of distinct pairs to size {@link #filter} for. */
    private final ObjectGoal<Long, P> sizeGoal;
    /** Deduplicates the pairs; held for the duration of the pass. */
    private ProbFilter filter;
    private Map<String, Counts> map;
    private Counts[] countsByPos;

    /**
     * Creates the goal, depending on the accession map, the loaded database and the sizing goal in
     * addition to the standard fasta-reader dependencies.
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
     * @param sizeGoal         the goal estimating the number of distinct pairs
     * @param deps             further goals this goal depends on
     */
    @SafeVarargs
    public DBQualityCountsGoal(P project, FTGoalKey key, ExecutionContext bundle, ObjectGoal<Set<RefSeqCategory>, P> categoriesGoal,
                               ObjectGoal<TaxNodeSelection, P> taxNodesGoal,
                               RefSeqFnaFilesDownloadGoal fnaFilesGoal,
                               ObjectGoal<Map<File, TaxTree.TaxIdNode>, P> additionalGoal,
                               ObjectGoal<AccessionMap, P> accessionMapGoal, ObjectGoal<Database, P> storeGoal,
                               ObjectGoal<Long, P> sizeGoal,
                               Goal<P>... deps) {
        super(project, key, bundle, categoriesGoal, taxNodesGoal, fnaFilesGoal, additionalGoal,
                accessionMapGoal, storeGoal, Goal.append(deps, sizeGoal));
        this.sizeGoal = sizeGoal;
    }

    /**
     * Re-reads the genomic fasta files, comparing their *k*-mers against the database to accumulate the
     * per-tax-id true-positive and positive counts, aggregates the counts up the selected ranks and
     * stores the resulting map as this goal's value.
     */
    @Override
    protected void doMakeThis() {
        try {
            prepare();
            long size = buildCounts();
            long bound = size;

            // How large the filter has to be. `size' above is the conservative bound: for every leaf it
            // sums the k-mers stored along its path to the root, i.e. it assumes each of them to occur
            // in every leaf below its node. That holds for a database whose k-mers sit close to the
            // leaves and fails badly for one whose weight is at a single high node -- `strepto' keeps
            // 222 million k-mers at the genus, on the path of every one of its several hundred leaves,
            // and the bound comes to some 10^11 entries. At ten bits each the filter alone would want
            // tens of gigabytes, and the run dies in newLargeGrid() before it reads a single base.
            //
            // The truth is roughly three times the stored k-mers, since a k-mer above the species is
            // hardly ever carried by more than three leaves. Rather than assume that factor, the pass
            // in `dbqualsize' sketches the pairs with HyperLogLog and counts them, exactly as `kmerindexsize'
            // does for the index filter -- same sketch class, same trade of one extra read for a
            // filter of the right size.
            GSConfigKey.BloomFilterSizing sizing =
                    (GSConfigKey.BloomFilterSizing) configValue(FTConfigKey.DB_QUALITY_FILTER_SIZING);
            if (sizing != GSConfigKey.BloomFilterSizing.UPPER_BOUND) {
                // Asking is what makes the sizing goal run: an ObjectGoal is a weak dependency, so
                // under the bound its pass over the sequences never happens at all.
                size = sizeGoal.get();
                if (getLogger().isInfoEnabled()) {
                    getLogger().info("Sizing the filter for " + size + " entries, the estimated number"
                            + " of distinct ones, against a bound of " + bound);
                }
            }

            // Using a blocked bloom filter here for more speed (identified the old Bloom filter as a
            // bottleneck). Allocated only now: under any sizing but the bound, `size' is what
            // `dbqualsize' counted, and that pass needs the store -- a k-mer that the database does not
            // hold forms no pair -- so the order is store, then size, then filter.
            filter = new BlockedBloomFilter(size);
            long bitSize = filter.getBitSize();
            if (getLogger().isInfoEnabled()) {
                getLogger().info("Filter size in MB: " + (bitSize / 8 / 1024 / 1024));
            }

            readFastas();

            long entries = 0;
            for (MyFastaReader reader : readers) {
                // The readers are done, so whatever the last (partial) batch still holds is looked up
                // and counted here, single-threaded, before their tallies are merged.
                reader.flushBatch();
                ((CountingReader) reader).mergeInto(countsByPos);
                entries += reader.entries;
            }
            if (getLogger().isInfoEnabled()) {
                getLogger().info("Filter entries: " + entries);
            }
            if (entries > 2 * size) {
                // Under the bound this means something went wrong. Under an estimate it may simply mean
                // the estimate fell short, and `auto' then says so with the remedy rather than leaving a
                // filter that was too small to have deduplicated reliably.
                if (sizing == GSConfigKey.BloomFilterSizing.UPPER_BOUND) {
                    if (getLogger().isErrorEnabled()) {
                        getLogger().error("Entries exceed filter size by over factor 2. Something went wrong!");
                    }
                } else {
                    throw new IllegalStateException("The sketched estimate of " + size + " distinct pairs fell short:"
                            + " " + entries + " were counted, over twice as many, so the filter was too small to"
                            + " deduplicate reliably. Set " + FTConfigKey.DB_QUALITY_FILTER_SIZING.getName() + "="
                            + GSConfigKey.BloomFilterSizing.UPPER_BOUND.getName() + " to use the conservative bound of "
                            + bound + " instead, if it can be allocated.");
                }
            }

            aggregateCounts(tree, map);
            set(map);
        } catch (IOException e) {
            throw new RuntimeException(e);
        } finally {
            map = null;
            countsByPos = null;
            filter = null;
            releaseAfterPass();
        }
    }

    /**
     * Builds the per-node tallies over the tree {@link #prepare()} loaded and returns the conservative
     * bound on how many pairs the pass will produce.
     * <p>
     * Separate from {@code prepare()} because only this goal has tallies: the sizing goal needs the
     * same tree, the same leaves and the same store, but counts nothing.
     *
     * @return the bound: for every leaf, the k-mers stored along its path to the root, summed
     */
    private long buildCounts() {
        Object2LongMap<String> stats = storeGoal.get().getStats();
        map = new HashMap<>();
        countsByPos = new Counts[nodeCount];
        // Estimate the filter size by summing up from species to root for each species in the DB.
        // It is a highly conservative estimate because k-mers on ranks above species are hardly
        // ever shared more than thrice (as found by measuring).
        long size = 0;
        for (SmallTaxTree.SmallTaxIdNode node : tree) {
            int pos = node.getPosition();
            boolean dataNode = leafByPos[pos];
            Counts counts = new Counts(dataNode, stats.getOrDefault(node.getTaxId(), 0L));
            if (dataNode) {
                // Count tp plus fp
                // Add k-mers from species upwards for each species:
                long pathSum = getPathSum(node, stats);
                counts.tpPlusFp = pathSum;
                size += pathSum;
            }
            map.put(node.getTaxId(), counts);
            countsByPos[pos] = counts;
        }
        return size;
    }

    @Override
    protected MyFastaReader newReader(AbstractRefSeqFastaReader.StringLong2DigitTrie contigsPerTaxid, AbstractRefSeqFastaReader.GenomeKeyTrie admittedGenomes) {
        return new CountingReader(contigsPerTaxid, admittedGenomes);
    }

    /**
     * Reader of this pass: it deduplicates every pair against the shared filter and tallies what
     * survives, per thread, into arrays that {@link #mergeInto} adds up afterwards.
     */
    protected class CountingReader extends MyFastaReader {
        /** Tallies of this reader alone, indexed by node position; merged by {@link #mergeInto}. */
        private final long[] tp;
        private final long[] tpPlusFn;
        private final long[] tpForNode;

        /**
         * Creates the reader and its tallies.
         *
         * @param contigsPerTaxid the trie counting contigs per tax id
         * @param admittedGenomes the shared set of genome keys admitted so far
         */
        CountingReader(AbstractRefSeqFastaReader.StringLong2DigitTrie contigsPerTaxid, AbstractRefSeqFastaReader.GenomeKeyTrie admittedGenomes) {
            super(contigsPerTaxid, admittedGenomes);
            tp = new long[nodeCount];
            tpPlusFn = new long[nodeCount];
            tpForNode = new long[nodeCount];
        }

        /**
         * Adds this reader's tallies to the shared per-node ones. Called once, after every reader has
         * finished, so nothing here needs a lock.
         *
         * @param target the per-position tallies to add to
         */
        void mergeInto(Counts[] target) {
            for (int pos = 0; pos < nodeCount; pos++) {
                if (tp[pos] != 0 || tpPlusFn[pos] != 0 || tpForNode[pos] != 0) {
                    target[pos].addFromReader(tp[pos], tpPlusFn[pos], tpForNode[pos]);
                }
            }
        }

        /**
         * Counts one (k-mer, leaf) pair unless it was seen before.
         *
         * @param kmer       the k-mer
         * @param leafPos    the position of the leaf it was read in
         * @param storedNode the node the database stores it at
         * @return whether the pair was new, i.e. whether it was counted
         */
        @Override
        protected boolean count(long kmer, int leafPos, SmallTaxTree.SmallTaxIdNode storedNode) {
            // Checks whether it's a duplicate under that leaf.
            if (!filter.putLong(KMerIndexFilterHelper.combine(kmer, leafPos))) {
                return false;
            }
            entries++;
            if (!leafByPos[leafPos]) {
                // The reader resolved a record to a node that isLeafNode() did not accept. Naming it
                // matters: the two halves disagree only when the chain of artificial nodes is not what
                // the flags imply, and the node's rank says which end is wrong.
                SmallTaxTree.SmallTaxIdNode leaf = nodeByPos[leafPos];
                throw new IllegalStateException("Must be a count for a leaf node here, but "
                        + leaf.getTaxId() + " (" + leaf.getName() + ") has rank " + leaf.getRank()
                        + " and is not the deepest artificial node on its branch.");
            }
            tpPlusFn[leafPos]++;
            // Is the stored node on the path from the leaf? The tree answers that from the two nodes'
            // depths: a node that is deeper than the leaf is ruled out without touching memory, and
            // otherwise exactly the depth difference is walked. Walking to the root instead - which is
            // what happens on every false negative, and on a refined tree of a few thousand levels
            // that is a few thousand steps - is what this replaces; the same change was the entire
            // measurable win in a profile of the matcher.
            if (tree.isAncestorOf(nodeByPos[leafPos], storedNode)) {
                // Stored node on the path from the leaf: true positive.
                tp[leafPos]++;
                tpForNode[storedNode.getPosition()]++;
            }
            return true;
        }
    }

    /**
     * Aggregates the per-node tallies up the taxonomy tree in two passes.
     * <p>
     * The first pass runs over the data taxa and adds each of them to all of its ancestors. It yields
     * the counts that refer to the path from a data taxon up to the <em>root</em>, which back the
     * per-taxon (unweighted) averages and the recalls.
     * <p>
     * The second pass covers the weighted path precision, which is defined relative to the node it is
     * reported for: for a node {@code m} it pools the counts of the paths from the data taxa under
     * {@code m} up to {@code m} itself, not up to the root. Those paths are exactly the nodes of the
     * subtree rooted at {@code m}, and since every k-mer of the database sits at exactly one node and
     * lies on the path of precisely the data taxa subordinate to it, pooling per node over the subtree
     * gives the same counts as pooling per data taxon over its path up to {@code m}. Reusing the first
     * pass instead would add the stretch from {@code m} to the root to every data taxon under {@code m}
     * and hence report a quantity that is not relative to {@code m} at all.
     *
     * @param nodes the nodes of the taxonomy tree, in any order
     * @param map   the tallies keyed by tax id, holding one entry for every node in {@code nodes}
     */
    protected static void aggregateCounts(Iterable<SmallTaxTree.SmallTaxIdNode> nodes, Map<String, Counts> map) {
        for (SmallTaxTree.SmallTaxIdNode node : nodes) {
            Counts counts = map.get(node.getTaxId());
            if (counts.isForLeaf()) {
                for (SmallTaxTree.SmallTaxIdNode ancestor = node; ancestor != null; ancestor = ancestor.getParent()) {
                    Counts c = map.get(ancestor.getTaxId());
                    // Leads to a weighted average:
                    // Node weight is proportional to positives of each aggregated node
                    c.aggregate(counts);
                }
            }
        }
        for (SmallTaxTree.SmallTaxIdNode node : nodes) {
            Counts counts = map.get(node.getTaxId());
            // Needs the leaf counts of the first pass, hence the separate loop.
            if (counts.getLeaves() > 0) {
                for (SmallTaxTree.SmallTaxIdNode ancestor = node; ancestor != null; ancestor = ancestor.getParent()) {
                    map.get(ancestor.getTaxId()).aggregateSubtree(counts);
                }
            }
        }
    }

    /**
     * Per-tax-id tally of true positives, true-positives-plus-false-positives and true-positives-plus-
     * false-negatives, plus aggregated precision/recall sums used to compute weighted and unweighted
     * averages.
     * <p>
     * Two frames of reference are kept apart here. The {@code tp} / {@code tpPlusFp} / {@code tpPlusFn}
     * counts refer to the path from a data taxon up to the <em>root</em> and, aggregated over the data
     * taxa under a node, back the unweighted averages and the recalls. The {@code subtree*} counts refer
     * to the paths from the data taxa under a node up to <em>that node</em> and back the weighted path
     * precision, which is defined relative to the subtree it is reported for.
     */
    public static class Counts implements Serializable {
        private static final long serialVersionUID = 1L;

        /**
         * Whether this tally belongs to a leaf, i.e. to the deepest artificial node on its branch --
         * see {@link DBQualityCountsGoal#isLeafNode}. Only leaf tallies are filled from the genomic
         * files; the others are aggregated from their leaves.
         */
        private final boolean forLeaf;

        /**
         * True positives plus false positives (all k-mers stored under this tax id).
         */
        private long tpPlusFp;
        /**
         * True positives (k-mers correctly stored under this tax id).
         */
        private long tp;
        /**
         * True positives plus false negatives (k-mers stored in the database that were also found in this tax id's genomes).
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
         * True positives summed over the subtree rooted at this node, i.e. over the paths from the data
         * taxa under this node up to this node itself. It backs the weighted path precision, which is
         * defined relative to this node rather than to the root.
         */
        private long subtreeTp;
        /**
         * True positives plus false positives summed over the subtree rooted at this node, i.e. over the
         * paths from the data taxa under this node up to this node itself.
         */
        private long subtreeTpPlusFp;
        /**
         * Sum of the per-k-mer precisions p(a) over the k-mers residing in the subtree rooted at this
         * node for which p is defined. It is the numerator of {@link #getSubtreePrecision()}.
         */
        private double subtreeKmerPrecisionSum;
        /**
         * Number of k-mers stored in the subtree rooted at this node, counting only those for which the
         * per-k-mer precision is defined. It is the denominator of {@link #getSubtreePrecision()}.
         */
        private long subtreeKmerSum;
        /**
         * Sum of the per-k-mer precisions p(a) over the k-mers of this subtree that are stored
         * <em>above</em> the data taxa. It is the numerator of {@link #getRestrictedSubtreePrecision()}.
         */
        private double subtreeKmerPrecisionSumAboveData;
        /**
         * Number of k-mers of this subtree stored above the data taxa, counting only those for which
         * the per-k-mer precision is defined. It is the denominator of
         * {@link #getRestrictedSubtreePrecision()}.
         */
        private long subtreeKmerSumAboveData;
        /**
         * True positives attributed to this node itself, i.e. k-mers read from a genome below it
         * that are stored in the database under exactly this tax id.
         */
        private long tpForNodePrecision;
        /**
         * Number of k-mers the database stores under exactly this tax id (not counting descendants).
         */
        private long kmerSumForNode;
        /**
         * Number of leaves aggregated into this tally, including this node itself if it is a leaf.
         */
        private int leaves;

        /**
         * Creates an empty tally with all counts set to zero.
         *
         * @param forLeaf        whether the tally belongs to a leaf node, i.e. the deepest artificial
         *                       node on its branch (see {@link DBQualityCountsGoal#isLeafNode})
         * @param kmerSumForNode the number of k-mers the database stores under exactly this tax id
         */
        public Counts(boolean forLeaf, long kmerSumForNode) {
            this.forLeaf = forLeaf;
            this.kmerSumForNode = kmerSumForNode;
        }

        /**
         * Returns whether this tally belongs to a leaf node.
         *
         * @return whether this tally belongs to a leaf node, i.e. the deepest artificial node on its
         * branch (see {@link DBQualityCountsGoal#isLeafNode})
         */
        public boolean isForLeaf() {
            return forLeaf;
        }

        /**
         * Returns the true positives attributed to this node itself.
         *
         * @return the number of read k-mers stored in the database under exactly this tax id
         */
        public long getTpForNodePrecision() {
            return tpForNodePrecision;
        }

        /**
         * Records one more k-mer of this node that was read from the genome of a data taxon underneath,
         * counted once per (k-mer, data taxon) pair.
         */
        void incTpForNodePrecision() {
            tpForNodePrecision++;
        }

        /**
         * Adds one reader's tallies for this node to it.
         * <p>
         * The readers count into arrays of their own and are merged here once they have all finished,
         * rather than incrementing these fields under a lock as they go. Two locks per
         * (k-mer, data taxon) pair is what that cost, and the pairs are counted in the billions on the
         * nodes every thread touches - the species a database is built for, and the data node beneath
         * it - so the contention fell on exactly the nodes that are hottest.
         *
         * @param tp                 true positives to add
         * @param tpPlusFn           true positives plus false negatives to add
         * @param tpForNodePrecision node-precision true positives to add
         */
        void addFromReader(long tp, long tpPlusFn, long tpForNodePrecision) {
            this.tp += tp;
            this.tpPlusFn += tpPlusFn;
            this.tpForNodePrecision += tpForNodePrecision;
        }

        /**
         * Returns the number of leaves aggregated into this tally.
         *
         * @return the number of leaves aggregated into this tally, including this node itself if it
         *         is a leaf
         */
        public int getLeaves() {
            return leaves;
        }

        /**
         * Returns the number of k-mers stored under this node alone.
         *
         * @return the number of k-mers the database stores under exactly this tax id
         */
        public long getKmerSumForNode() {
            return kmerSumForNode;
        }

        /**
         * Returns the node precision of this node, i.e. the mean of the per-k-mer precisions p(a) over
         * the k-mers stored under this node alone. The precision of a single k-mer a is the fraction
         * c(a) / |D| of the data taxa under a's node whose genomes actually contain a, so averaging it
         * over this node's k-mers gives {@code tpForNodePrecision / (leaves * kmerSumForNode)}.
         * <p>
         * The value is defined only for nodes that hold k-mers and have at least one data taxon
         * underneath; for all others it is {@link Double#NaN}. It is deliberately not extended to
         * those by a convention such as {@code 1}: a node storing no k-mer makes no statement about
         * k-mer placement, and a perfect score would let it raise {@link #getSubtreePrecision()}.
         * For the same reason no Laplace correction is applied - a virtual k-mer per node would enter
         * numerator and denominator alike, raising a ratio that is at most one anyway, by an amount
         * growing with the number of nodes of the taxonomy.
         *
         * @return the precision of this node's own k-mers, or {@link Double#NaN} if the node holds no
         * k-mers or has no data taxon underneath
         */
        public double getNodePrecision() {
            if (leaves == 0 || kmerSumForNode == 0) {
                // Nothing to be precise about: no k-mers stored here, or no data-taxon descendants
                // counted (e.g. OTHER nodes). The measure is undefined rather than perfect.
                return Double.NaN;
            }
            return ((double) tpForNodePrecision) / (leaves * kmerSumForNode);
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
         * Returns the weighted average path precision of the subtree rooted at this node, i.e. the
         * pooled true positives over the pooled positives of the paths from the data taxa under this
         * node up to this node itself. Note that this is relative to this node and not to the root:
         * the k-mers stored above this node are shared by all data taxa underneath and say nothing
         * about how well the subtree itself is refined.
         * <p>
         * The value is deliberately not Laplace-corrected. A virtual k-mer per node would enter
         * numerator and denominator alike and hence raise a ratio that is at most one anyway, by an
         * amount that grows with the number of nodes of the subtree. As the refinement adds nodes to
         * the taxonomy, a corrected measure would reward a refined tree over an unrefined one for its
         * node count alone. Without the correction the value depends on the placement of the k-mers
         * only: it equals the sum of c(a) over the sum of |D| of the holding node, taken over the
         * k-mers of the subtree, and inserting an inner node changes neither of the two.
         *
         * @return the weighted path precision of this node's subtree, or {@link Double#NaN} if the
         * subtree holds no k-mer at all
         */
        public double getPrecision() {
            if (subtreeTpPlusFp == 0) {
                return Double.NaN;
            }
            return ((double) subtreeTp) / subtreeTpPlusFp;
        }

        /**
         * Returns the subtree precision restricted to the k-mers stored above the data taxa, i.e.
         * averaged over those k-mers alone. This is where a refinement can act: a k-mer at a data
         * taxon is fixed at a precision of one by construction.
         *
         * @return the restricted subtree precision, or {@link Double#NaN} if the subtree stores no
         * k-mer above its data taxa
         */
        public double getRestrictedSubtreePrecision() {
            if (subtreeKmerSumAboveData == 0) {
                return Double.NaN;
            }
            return subtreeKmerPrecisionSumAboveData / subtreeKmerSumAboveData;
        }

        /**
         * Returns the number of k-mers backing {@link #getRestrictedSubtreePrecision()}.
         *
         * @return the number of k-mers of this subtree stored above its data taxa
         */
        public long getSubtreeKMersAboveData() {
            return subtreeKmerSumAboveData;
        }

        /**
         * Returns the subtree precision of the subtree rooted at this node, i.e. the mean of the
         * per-k-mer precisions p(a) over the k-mers residing in that subtree. It is the same average
         * as {@link #getNodePrecision()}, just taken over the k-mers of a whole subtree instead of
         * those of a single node.
         * <p>
         * Equivalently, it is the mean of the node precisions within the subtree, each weighted by that
         * node's share of the subtree's k-mers; a node without k-mers has no share and simply does not
         * occur. Since every k-mer counts equally, inserting an inner node into the taxonomy - which
         * leaves every p(a) untouched - cannot change the value. Moving a k-mer downwards, in turn,
         * leaves c(a) untouched (path correctness fixes it) while shrinking |D| of its node, so p(a)
         * can only rise: the value responds to nothing but the actual relocation of k-mers, and it
         * responds to it in the right direction. Note that the k-mers residing above this node do not
         * enter, as they are shared by all data taxa underneath.
         *
         * @return the subtree precision of this node's subtree, or {@link Double#NaN} if the subtree
         * holds no k-mer for which the per-k-mer precision is defined
         */
        public double getSubtreePrecision() {
            if (subtreeKmerSum == 0) {
                return Double.NaN;
            }
            return subtreeKmerPrecisionSum / subtreeKmerSum;
        }

        /**
         * Returns the number of k-mers backing {@link #getSubtreePrecision()}, i.e. those residing in
         * the subtree rooted at this node for which the per-k-mer precision is defined.
         *
         * @return the number of k-mers the subtree precision averages over
         */
        public long getSubtreeKmerSum() {
            return subtreeKmerSum;
        }

        /**
         * Returns the true positives pooled over the subtree rooted at this node.
         *
         * @return the true positives of the paths from the data taxa under this node up to this node
         */
        public long getSubtreeTp() {
            return subtreeTp;
        }

        /**
         * Returns the true positives plus false positives pooled over the subtree rooted at this node.
         *
         * @return the positives of the paths from the data taxa under this node up to this node
         */
        public long getSubtreeTpPlusFp() {
            return subtreeTpPlusFp;
        }

        /**
         * Returns the path precision of a single data taxon, i.e. relative to the root rather than to a
         * subtree. It backs the unweighted averages.
         *
         * @return the precision {@code tp / (tp + fp)}
         */
        private double getRawPrecision() {
            return ((double) tp / tpPlusFp);
        }

        // Unweighted

        /**
         * Returns the unweighted average precision, i.e. the mean of the per-taxon path precisions of
         * the data taxa under this node, each taken over the path up to the root.
         *
         * @return the unweighted average precision over aggregated nodes, or this node's own path
         * precision if there was no aggregation, i.e. if this node is a data taxon itself
         */
        public double getAvgPrecision() {
            if (aggPrecisionSum == 0) {
                // No data taxon was aggregated, so this node is one itself and carries its own path
                // precision - which, unlike the weighted average, refers to the path up to the root.
                return getRawPrecision();
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
            leaves++;
            if (!isForLeaf()) {
                tp += counts.tp;
                tpPlusFp += counts.tpPlusFp;
                tpPlusFn += counts.tpPlusFn;
                aggregations++;
                aggPrecisionSum += counts.getAvgPrecision();
                aggRecallSum += counts.getAvgRecall();
            }
        }

        /**
         * Adds the contribution of a single node of this node's subtree - which may be this node itself
         * - to the pooled counts backing {@link #getPrecision()}. Every k-mer of the database sits at
         * exactly one node and lies on the path of precisely those data taxa that are subordinate to
         * that node, so summing per node over the subtree yields the same counts as summing per data
         * taxon over its path up to this node.
         *
         * @param counts the tally of a node of this node's subtree, which must have at least one data
         *               taxon underneath
         */
        private void aggregateSubtree(Counts counts) {
            // Nodes without k-mers contribute nothing to any of these sums, which is what makes both
            // averages insensitive to the number of nodes the taxonomy provides.
            subtreeTp += counts.tpForNodePrecision;
            subtreeTpPlusFp += counts.leaves * counts.kmerSumForNode;
            if (counts.kmerSumForNode > 0) {
                // The sum of p(a) = c(a) / |D_n| over the k-mers a stored at n. The caller guarantees
                // counts.leaves > 0, so p is defined for all of them.
                double precisionSum = ((double) counts.tpForNodePrecision) / counts.leaves;
                subtreeKmerPrecisionSum += precisionSum;
                subtreeKmerSum += counts.kmerSumForNode;
                if (!counts.forLeaf) {
                    // A leaf has no data taxa under it but itself after the tree transformation,
                    // so |D_n| = 1 and every k-mer stored there has p(a) = 1 whatever the database
                    // looks like. Since such
                    // k-mers are the bulk of a database, they dominate the average while being
                    // incapable of improvement; the restricted sums leave them out.
                    subtreeKmerPrecisionSumAboveData += precisionSum;
                    subtreeKmerSumAboveData += counts.kmerSumForNode;
                }
            }
        }
    }
}
