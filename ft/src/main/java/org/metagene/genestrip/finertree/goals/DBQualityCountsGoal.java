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
import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.probfilter.BlockedBloomFilter;
import org.metagene.genestrip.probfilter.ProbFilter;
import org.metagene.genestrip.finertree.FTGoalKey;
import org.metagene.genestrip.finertree.FTProject;
import org.metagene.genestrip.finertree.FTConfigKey;
import org.metagene.genestrip.finertree.probfilter.DistinctPairSketch;
import org.metagene.genestrip.finertree.probfilter.KMerIndexFilterHelper;
import org.metagene.genestrip.finertree.refseq.AbstractUpdateFastaReader;
import org.metagene.genestrip.genbank.AssemblySummaryReader;
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
import org.metagene.genestrip.store.RadixKMerStore;
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
    /**
     * Number of k-mers a reader buffers before looking them up in one batch. Sized to keep enough
     * independent cache-missing lookups in flight for the memory-level parallelism to pay, without the
     * per-batch bookkeeping outweighing it - the same value {@code AbstractKMerIndexGoal} uses.
     */
    protected static final int BATCH_SIZE = 128;

    private final ObjectGoal<AccessionMap, P> accessionMapGoal;
    private final ObjectGoal<Database, P> storeGoal;

    private SmallTaxTree tree;
    private KMerStore<SmallTaxTree.SmallTaxIdNode> kMerSortedArray;
    private ProbFilter filter;
    /**
     * Non-null only during the sizing pass, when the readers sketch their (k-mer, leaf) pairs instead
     * of counting them. See the sizing block of {@link #doMakeThis()}.
     */
    private DistinctPairSketch sketch;
    private Map<String, Counts> map;
    private List<MyFastaReader> readersList;
    // The reading path addresses a node by its dense position and never by its tax id: the map above
    // is keyed by a String, and looking a Counts up in it twice per (k-mer, data taxon) pair hashed a
    // string a few billion times per run. SmallTaxTree numbers its nodes densely for exactly this
    // (see SmallTaxTree#getNodeCount), so an array does the same job with an indexed read.
    private SmallTaxTree.SmallTaxIdNode[] nodeByPos;
    private Counts[] countsByPos;
    private boolean[] leafByPos;
    private int nodeCount;

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
        GSProject project = getProject();
        // Data nodes, and not merely one of the three kinds of artificial node. What a leaf is has to
        // agree between the fill and this measure, and `dataNodes' is what makes that agreement
        // simple: with it on, ReworkingStoreFastaReader.reworkNode() files every genome at a DATA
        // node or deeper, so no taxonomy node ever holds a genome's k-mers and a node without
        // children is exactly a node the fill filed a genome at. isLeafNode() is then one test and
        // needs to know nothing about ranks.
        //
        // Requiring it is not what used to shut this goal out of a database refined below the
        // species -- that was the second half of the old guard, which REJECTED `fileNodes' while
        // only a DATA node could be recognised as a leaf. Both halves went at once; only the first
        // is coming back. `fileNodes' and `idNodes' stay free, and the `cdiff' project of
        // ft-db-exp2, which needs file nodes because the taxonomy supplies no children below the
        // species, has data nodes on as every project here does.
        if (!project.booleanConfigValue(GSConfigKey.DATA_NODES)) {
            throw new IllegalStateException("This goal requires data nodes (dataNodes=true)");
        }

        try {
            map = new HashMap<>();
            tree = storeGoal.get().getTaxTree();
            Object2LongMap<String> stats = storeGoal.get().getStats();
            // Estimate the filter size by summing up from species to root for each species in the DB.
            // It is a highly conservative estimate because k-mers on ranks above species are hardly
            // ever shared more than thrice (as found by measuring).
            nodeCount = tree.getNodeCount();
            nodeByPos = new SmallTaxTree.SmallTaxIdNode[nodeCount];
            countsByPos = new Counts[nodeCount];
            leafByPos = new boolean[nodeCount];
            long size = 0;
            for (SmallTaxTree.SmallTaxIdNode node : tree) {
                boolean dataNode = isLeafNode(node);
                Counts counts = new Counts(dataNode, stats.getOrDefault(node.getTaxId(), 0L));
                if (dataNode) {
                    // Count tp plus fp
                    // Add k-mers from species upwards for each species:
                    long pathSum = getPathSum(node, stats);
                    counts.tpPlusFp = pathSum;
                    size += pathSum;
                }
                map.put(node.getTaxId(), counts);
                int pos = node.getPosition();
                nodeByPos[pos] = node;
                countsByPos[pos] = counts;
                leafByPos[pos] = dataNode;
            }
            kMerSortedArray = storeGoal.get().convertKMerStore();

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
            // below sketches the pairs with HyperLogLog and counts them, exactly as `kmerindexsize'
            // does for the index filter -- same sketch class, same trade of one extra read for a
            // filter of the right size.
            GSConfigKey.BloomFilterSizing sizing =
                    (GSConfigKey.BloomFilterSizing) configValue(FTConfigKey.DB_QUALITY_FILTER_SIZING);
            long bound = size;
            if (sizing != GSConfigKey.BloomFilterSizing.UPPER_BOUND) {
                sketch = new DistinctPairSketch(1);
                readersList = new ArrayList<>();
                try {
                    readFastas();
                    for (MyFastaReader reader : readersList) {
                        reader.flushBatch();
                    }
                    size = sketch.estimate();
                } finally {
                    sketch = null;
                    readersList = null;
                }
                if (getLogger().isInfoEnabled()) {
                    getLogger().info("Estimated distinct filter entries: " + size
                            + ", against a bound of " + bound);
                }
            }

            // Using a blocked bloom filter here for more speed (identified the old Bloom filter as a
            // bottleneck). Allocated only now: under any sizing but the bound, `size' is what the pass
            // above counted, and that pass needs the store -- a k-mer that the database does not hold
            // forms no pair -- so the order is store, then size, then filter.
            filter = new BlockedBloomFilter(size);
            long bitSize = filter.getBitSize();
            if (getLogger().isInfoEnabled()) {
                getLogger().info("Filter size in MB: " + (bitSize / 8 / 1024 / 1024));
            }

            readersList = new ArrayList<>();
            readFastas();

            long entries = 0;
            for (MyFastaReader reader : readersList) {
                // The readers are done, so whatever the last (partial) batch still holds is looked up
                // and counted here, single-threaded, before their tallies are merged.
                reader.flushBatch();
                reader.mergeInto(countsByPos);
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
            tree = null;
            kMerSortedArray = null;
            filter = null;
            sketch = null;
            readersList = null;
            nodeByPos = null;
            countsByPos = null;
            leafByPos = null;
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
        MyFastaReader reader = new MyFastaReader(intConfigValue(GSConfigKey.FASTA_LINE_SIZE_BYTES),
                taxNodesGoal.get(),
                isIncludeRefSeqFna() ? accessionMapGoal.get() : null,
                intConfigValue(GSConfigKey.KMER_SIZE),
                intConfigValue(GSConfigKey.MAX_GENOMES_PER_TAXID),
                (Rank) configValue(GSConfigKey.MAX_GENOMES_PER_TAXID_RANK),
                longConfigValue(GSConfigKey.MAX_KMERS_PER_TAXID),
                intConfigValue(GSConfigKey.MAX_DUST),
                intConfigValue(GSConfigKey.KMER_SAMPLING),
                booleanConfigValue(GSConfigKey.ASSEMBLY_ACCESSIONS_ONLY),
                regionsPerTaxid,
                booleanConfigValue(GSConfigKey.ENABLE_LOWERCASE_BASES));
        readersList.add(reader);
        return reader;
    }

    /**
     * Fasta reader that, for each *k*-mer read from a genome, checks whether it is stored in the
     * database and whether the stored node lies on the path from the read's leaf node, updating the
     * per-tax-id counts accordingly while deduplicating via the bloom filter.
     */
    protected class MyFastaReader extends AbstractUpdateFastaReader
            implements RadixKMerStore.BatchValueConsumer<SmallTaxTree.SmallTaxIdNode> {
        /** Number of distinct (k-mer, leaf node) pairs this reader added to the dedup filter. */
        protected long entries;

        /**
         * Buffers for the batched store lookup, or {@code null} when it cannot be used. The store
         * lookup is memory-latency bound, and a batch lets many of its cache misses overlap - see
         * {@link RadixKMerStore#getBatch}.
         */
        private final RadixKMerStore.BatchBuffers batch;
        /** Tallies of this reader alone, indexed by node position; merged by {@link #mergeInto}. */
        private final long[] tp;
        private final long[] tpPlusFn;
        private final long[] tpForNode;
        /** The current region's leaf and its position, resolved once per region rather than per k-mer. */
        private SmallTaxTree.SmallTaxIdNode cachedLeaf;
        private int cachedLeafPos = -1;

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
         * @param kMerSampling               the step size between stored k-mers
         * @param assemblyAccessionsOnly    whether only complete genomes are considered
         * @param regionsPerTaxid        the trie counting regions per tax id
         * @param enableLowerCaseBases   whether lower-case bases are treated as regular bases
         */
        public MyFastaReader(int bufferSize, Set<TaxTree.TaxIdNode> taxNodes, AccessionMap accessionMap,
                             int k, int maxGenomesPerTaxId, Rank maxGenomesPerTaxIdRank, long maxKmersPerTaxId, int maxDust, int kMerSampling, boolean assemblyAccessionsOnly, StringLong2DigitTrie regionsPerTaxid, boolean enableLowerCaseBases) {
            super(bufferSize, taxNodes, accessionMap, k, maxGenomesPerTaxId, maxGenomesPerTaxIdRank, maxKmersPerTaxId, maxDust, kMerSampling, assemblyAccessionsOnly, regionsPerTaxid, enableLowerCaseBases, booleanConfigValue(GSConfigKey.ID_NODES), booleanConfigValue(GSConfigKey.FILE_NODES), booleanConfigValue(GSConfigKey.DATA_NODES));
            entries = 0;
            tp = new long[nodeCount];
            tpPlusFn = new long[nodeCount];
            tpForNode = new long[nodeCount];
            // Batched only while no per-taxon limit binds. A batched k-mer is counted after
            // handleStore() has already returned, so the return value can no longer say whether it
            // was, and that value feeds kmersInRegion - which endRegion() adds to the per-taxon
            // counters that maxGenomesPerTaxid and maxKMersPerTaxid are enforced from. At the
            // defaults neither binds and nothing reads those counters back; with either set, the
            // one-at-a-time path keeps the accounting exact.
            boolean unlimited = maxGenomesPerTaxId == Integer.MAX_VALUE && maxKmersPerTaxId == Long.MAX_VALUE;
            batch = (kMerSortedArray instanceof RadixKMerStore && unlimited)
                    ? new RadixKMerStore.BatchBuffers(BATCH_SIZE) : null;
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
         * Looks the buffered k-mers up in one batch and counts those the database holds. Called when
         * the buffer fills and, for the trailing ones, after all fastas have been read.
         */
        protected void flushBatch() {
            if (batch != null && !batch.isEmpty()) {
                ((RadixKMerStore<SmallTaxTree.SmallTaxIdNode>) kMerSortedArray).getBatch(batch, this);
            }
        }

        /**
         * Resolves the region's leaf as usual and caches its position, so that the reading path costs
         * a field read per k-mer instead of a lookup. The leaf changes per region; the k-mers of a
         * region are legion.
         */
        @Override
        protected void updateLeafNode() {
            super.updateLeafNode();
            if (leafNode != cachedLeaf) {
                cachedLeaf = leafNode;
                cachedLeafPos = -1;
                if (leafNode != null) {
                    if (leafNode.storeIndex < 0) {
                        // A leaf the store knows no value for. Database.initStoreIndex() gives every
                        // node of the tree one, so this says the tree and the store have come apart -
                        // worth a word, and once per region rather than once per k-mer.
                        if (getLogger().isWarnEnabled()) {
                            getLogger().warn("No kmer-index for taxid " + leafNode.getTaxId() + " found.");
                        }
                    } else {
                        cachedLeafPos = leafNode.getPosition();
                    }
                }
            }
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
        protected boolean handleStore(long kmer) {
            if (cachedLeafPos < 0) {
                return false;
            }
            if (batch != null) {
                // The leaf is only known while reading, so its position travels with the k-mer as the
                // batch's payload and the counting happens in accept() once the whole batch has been
                // looked up. Whether this k-mer will be counted is not known yet, hence false.
                if (batch.add(kmer, cachedLeafPos)) {
                    flushBatch();
                }
                return false;
            }
            SmallTaxTree.SmallTaxIdNode storedNode = kMerSortedArray.getLong(kmer, null);
            // There may be no corresponding node in the database:
            if (storedNode == null) {
                return false;
            }
            return count(kmer, cachedLeafPos, storedNode);
        }

        /**
         * Counts one k-mer of a flushed batch, which by construction the database holds - so the
         * {@code null} check of the unbatched path is implicit here.
         *
         * @param kmer       the k-mer that was looked up
         * @param leafPos    the position of the leaf it was read in, as buffered with it
         * @param storedNode the node the database stores it at
         */
        @Override
        public void accept(long kmer, int leafPos, SmallTaxTree.SmallTaxIdNode storedNode) {
            count(kmer, leafPos, storedNode);
        }

        /**
         * Counts one (k-mer, leaf) pair unless it was seen before.
         *
         * @param kmer       the k-mer
         * @param leafPos    the position of the leaf it was read in
         * @param storedNode the node the database stores it at
         * @return whether the pair was new, i.e. whether it was counted
         */
        private boolean count(long kmer, int leafPos, SmallTaxTree.SmallTaxIdNode storedNode) {
            if (sketch != null) {
                // The sizing pass: how many distinct pairs there are is the whole question, so the pair
                // is sketched and nothing is counted. The filter does not exist yet -- its size is what
                // this pass is for.
                sketch.record(KMerIndexFilterHelper.combine(kmer, leafPos));
                return true;
            }
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
     * Whether the given node is where a genomic file's k-mers come to rest, and therefore the unit
     * this goal's measures are taken over.
     * <p>
     * The database fill nests the artificial nodes: {@link org.metagene.genestrip.refseq.ReworkingStoreFastaReader#reworkNode()}
     * descends a tax id into its {@link Rank#DATA} child, that into a {@link Rank#FILE} child, and
     * that into a {@link Rank#ID} child, as far as {@code dataNodes}, {@code fileNodes} and
     * {@code idNodes} are enabled, and stores the k-mers at whichever it ends on. Reading the fastas
     * back, {@link AbstractUpdateFastaReader#updateLeafNode()} walks the same chain from the other
     * end -- ID, then FILE, then DATA, returning at the first that exists. Both therefore land on the
     * <em>deepest</em> of the three, which is what this method identifies: an origin-rank node with
     * no origin-rank child.
     * <p>
     * Testing for {@link Rank#DATA} alone, as this did before, is only equivalent while
     * {@code fileNodes} and {@code idNodes} are both off. With either on, the data node becomes an
     * empty intermediate holding no k-mers of its own while the reader resolves records to the file
     * or id node below it, and the two halves of this goal disagree about what a leaf is -- which
     * {@link DBQualityCountsGoal.MyFastaReader#handleStore(long)} catches and turns into an {@link IllegalStateException}.
     * <p>
     * {@link Rank#REFINED} is deliberately not an origin rank. A refined node is inserted by the
     * refinement <em>above</em> the origin nodes and holds the k-mers it moved down there, so it is
     * internal in exactly the way a taxonomy node is, and the measures restricted to what sits above
     * the data (see {@link Counts#aggregateSubtree}) have to keep counting it.
     *
     * @param node the node to test
     * @return whether the node is the deepest artificial node on its branch
     */
    protected static boolean isLeafNode(SmallTaxTree.SmallTaxIdNode node) {
        // A node the fill filed a genome at, which with data nodes on (see doMakeThis) is exactly a
        // node without children: every genome gets a DATA node or something below it, and whatever
        // the refinement inserts between a node and its original children gives that node children
        // and so keeps it internal. Asking about the children's ranks instead is what got this
        // wrong: after a refinement a data node's file nodes are no longer its children, the data
        // node passed for a leaf, and its k-mers -- the ones a refinement has the most to gain on --
        // dropped out of every average restricted to what sits above the data taxa. On cdiff that
        // alone lifted the reported sp* from 0.236 to 0.338 with no k-mer moving.
        SmallTaxTree.SmallTaxIdNode[] subNodes = node.getSubNodes();
        if (subNodes != null && subNodes.length != 0) {
            return false;
        }
        // ... with one exception: the "OTHER" placeholder the refinement inserts is childless but is
        // not a data taxon. Section "Data taxa and path correctness" defines one as a taxon with a
        // complete genome directly associated whose k-mers are stored in the database, whereas OTHER
        // stands for exactly the taxa that are *not* in the database and merely caused a k-mer to be
        // pushed above the species during the LCA update. It therefore never receives a k-mer -- in
        // the six databases of the paper all 1,876 of them are empty -- and counting it as a leaf
        // added one to |D_n| for every node the refinement touched, dividing p(a) = c(a)/|D_nu(a)|
        // accordingly without a single k-mer having moved. Where a genus held a single species that
        // was a halving: vineyard's Coniella, Pseudopezicula and Trichothecium each reported a
        // restricted subtree precision of exactly 0.5 against 1.0 before the refinement.
        //
        // The test is structural rather than by name. UpdateStoreGoal.createNode gives the REFINED
        // rank to two kinds of node: internal dendrogram nodes, which always have two children, and
        // the OTHER bucket, which is a leaf. A childless REFINED node is therefore the placeholder
        // and nothing else -- checked against all six databases, where the two sets coincide exactly.
        // getRankOrdinal() rather than getRank(), which is null for a rank the Rank enum does not
        // know; REFINED always has one, but the null-safe accessor keeps the guard honest.
        return node.getRankOrdinal() != Rank.REFINED.ordinal();
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
         * Returns this subtree's precision: the k-mer-weighted mean of the per-node precisions over
         * all of the subtree's k-mers, or {@link Double#NaN} if the subtree holds none.
         *
         * @return the subtree precision, or {@code NaN} if the subtree holds no k-mers
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
