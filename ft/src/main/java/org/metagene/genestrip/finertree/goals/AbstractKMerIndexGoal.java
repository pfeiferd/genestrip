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

import org.metagene.genestrip.ExecutionContext;
import org.metagene.genestrip.GSConfigKey;
import org.metagene.genestrip.finertree.FTConfigKey;
import org.metagene.genestrip.finertree.FTProject;
import org.metagene.genestrip.finertree.probfilter.KMerIndexFilterHelper;
import org.metagene.genestrip.finertree.refseq.AbstractUpdateFastaReader;
import org.metagene.genestrip.goals.refseq.FastaReaderGoal;
import org.metagene.genestrip.goals.refseq.RefSeqFnaFilesDownloadGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.GoalKey;
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
import org.metagene.genestrip.tax.TaxIdCollector;
import org.metagene.genestrip.tax.TaxTree;

import java.io.File;
import java.util.*;

/**
 * What the two goals walking the k-mer index have in common: both work out which tax nodes are
 * selected for refinement, both read every reference sequence, and both turn each k-mer stored below
 * such a node into a (k-mer, leaf index) pair. They part company only over what becomes of that
 * pair - one counts the pairs to size a filter, the other fills the filter - which is what
 * {@link #record(long)} decides.
 *
 * @param <T> the type of result the concrete goal produces
 * @param <P> the concrete {@link FTProject} type this goal operates on
 */
public abstract class AbstractKMerIndexGoal<T, P extends FTProject> extends FastaReaderGoal<T, P>
        implements Goal.LogHeapInfo {
    /**
     * Marker index used for k-mers that belong to no relevant leaf tax node, i.e. that must be
     * attributed to the "OTHER" bucket rather than to a concrete leaf.
     */
    public static final int OTHER_VALUE = Integer.MAX_VALUE;

    /**
     * Number of k-mers a {@link AbstractKMerIndexGoal.MyFastaReader} buffers before looking them up in one batch via
     * {@link RadixKMerStore#getBatch(RadixKMerStore.BatchBuffers, RadixKMerStore.BatchValueConsumer)}.
     * The store lookup is memory-latency bound, so batching lets many of its cache misses overlap.
     */
    protected static final int BATCH_SIZE = 128;

    private final ObjectGoal<AccessionMap, P> accessionMapGoal;
    private final ObjectGoal<Database, P> storeGoal;
    private final ObjectGoal<TaxTree, P> taxTreeGoal;

    /** The refinement positions, i.e. the ranks and tax ids at which the tree is refined. */
    protected final List<FTConfigKey.RefinementPosition> refinementIntervals;

    /** The k-mer store of the database being indexed, held while this goal is being made. */
    protected KMerStore<SmallTaxTree.SmallTaxIdNode> kmerStore;
    private SmallTaxTree smallTaxTree;
    private Set<TaxTree.TaxIdNode> relevantNodes;
    /**
     * Bit set marking the stored tax nodes selected for refinement, precomputed once so that
     * {@link MyFastaReader#handleStore(long)} avoids re-evaluating {@link FTConfigKey.RefinementPosition}
     * for every k-mer. A node is addressed by its store index, so the lookup is a plain bit test.
     */
    private BitSet refinementBits;

    /** The readers of the current pass, one per thread. */
    protected final List<MyFastaReader> readers = new ArrayList<>();

    /**
     * Creates the goal.
     *
     * @param project          the project this goal belongs to
     * @param goalKey          the key identifying the concrete goal
     * @param bundle           the execution context providing the number of worker threads
     * @param categoriesGoal   goal supplying the RefSeq categories to process
     * @param taxNodesGoal     goal supplying the tax nodes to be considered
     * @param taxTreeGoal      goal supplying the (large) tax tree
     * @param fnaFilesGoal     goal supplying the downloaded RefSeq fna files
     * @param additionalGoal   goal supplying additional fasta files mapped to their tax node
     * @param accessionMapGoal goal supplying the accession-to-tax-node mapping
     * @param storeGoal        goal supplying the database (and its k-mer store) to index
     * @param deps             further goals this goal depends on
     */
    @SafeVarargs
    public AbstractKMerIndexGoal(P project, GoalKey goalKey, ExecutionContext bundle,
                                 ObjectGoal<Set<RefSeqCategory>, P> categoriesGoal,
                                 ObjectGoal<Set<TaxTree.TaxIdNode>, P> taxNodesGoal,
                                 ObjectGoal<TaxTree, P> taxTreeGoal, RefSeqFnaFilesDownloadGoal fnaFilesGoal,
                                 ObjectGoal<Map<File, TaxTree.TaxIdNode>, P> additionalGoal,
                                 ObjectGoal<AccessionMap, P> accessionMapGoal, ObjectGoal<Database, P> storeGoal,
                                 Goal<P>... deps) {
        // true, and not the configured `refseq.filldb': this pass reads the RefSeq release whether or
        // not the database was filled from it, for the same reason the lowest-common-ancestor update
        // does. What keeps a k-mer from being pushed down to a child is a genome that carries it and
        // sits elsewhere - and when the release did not fill the database, those genomes are exactly
        // the ones only the release has. Skipping it would leave them unseen, no OTHER slot would be
        // set for them, and every such k-mer would be attributed to whichever children happen to be in
        // the database - a specificity the genomes outside it contradict.
        super(project, goalKey, bundle, categoriesGoal, taxNodesGoal, fnaFilesGoal, additionalGoal, true,
                Goal.append(deps, taxTreeGoal, accessionMapGoal, storeGoal));
        this.storeGoal = storeGoal;
        this.accessionMapGoal = accessionMapGoal;
        this.taxTreeGoal = taxTreeGoal;
        refinementIntervals = (List<FTConfigKey.RefinementPosition>) configValue(FTConfigKey.REFINEMENT_POSITIONS);
    }

    /**
     * Does what both goals must do before they can read anything: determines the tax nodes selected
     * for refinement and loads the k-mer store they are looked up in.
     */
    protected void prepare() {
        TaxIdCollector collector = new TaxIdCollector(taxTreeGoal.get());
        // Quite inefficient but should be good enough at this place.
        Set<TaxTree.TaxIdNode> nodesWithRank = new HashSet<>();
        for (TaxTree.TaxIdNode node : taxNodesGoal.get()) {
            while (node != null) {
                if (FTConfigKey.RefinementPosition.getMatchingNodeFor(node, refinementIntervals) != null) {
                    nodesWithRank.add(node);
                }
                node = node.getParent();
            }
        }
        // Include subnodes from all nodes where we have ranks to refine:
        relevantNodes = collector.withDescendants(nodesWithRank, (Rank) configValue(GSConfigKey.RANK_COMPLETION_DEPTH));
        smallTaxTree = storeGoal.get().getTaxTree();
        kmerStore = storeGoal.get().convertKMerStore();
        // Precompute which stored nodes are selected for refinement, so the per-k-mer hot path in
        // handleStore() only needs a bit test.
        initRefinementBits(kmerStore.getNKmersPerTaxid().keySet());
    }

    /**
     * Releases what {@link #prepare()} loaded. The reader threads are not this method's business:
     * a pass ends its own consumers when it has read everything.
     */
    protected void releaseAfterPass() {
        relevantNodes = null;
        smallTaxTree = null;
        kmerStore = null;
        refinementBits = null;
        readers.clear();
    }

    /**
     * Records one (k-mer, leaf index) pair, already combined into a single value.
     *
     * @param hash the combined pair
     * @return whether the pair was one this goal had not recorded before
     */
    protected abstract boolean record(long hash);

    /**
     * Returns the slot a (k-mer, leaf) pair is recorded under: the store index of the <em>direct
     * child</em> of the node being refined that the leaf lies under, or {@link #OTHER_VALUE} when
     * there is none.
     * <p>
     * Not the leaf's own index, deep in the tree, although that is where the pair comes from. The
     * reassignment asks the filter once per direct child subtree of the node being refined and once
     * for {@link #OTHER_VALUE}, and for nothing else - see
     * {@link KMerStoreWorkGoal}, whose {@code bits} array has exactly that many slots and is the
     * whole of what the clustering, the histogram and the store rewrite ever read. Which genome
     * below a child carries the k-mer is a distinction nothing acts on, so recording it costs
     * entries and buys nothing: a k-mer that all 3,500 assemblies of a species carry needs one entry
     * per child subtree, not one per assembly. Registering at the child is also what makes
     * {@code KMerIndexBloomGoal.conservativeEstimate()} a bound in fact and not just in intent - the
     * pairs a node can produce are then its children plus OTHER, which is precisely what that sum
     * counts.
     * <p>
     * {@link #OTHER_VALUE} covers two cases, and the walk up the tree answers both at once. The leaf
     * may <em>be</em> the node being refined - a genome read from a source that got no file node of
     * its own resolves to the taxon itself - and the pair then belongs in the OTHER slot, or it is
     * recorded where nobody looks and the k-mer would appear not to be carried by that genome at
     * all. The leaf may also lie outside the node's subtree altogether, which cannot happen while
     * the stored node is the lowest common ancestor of every region carrying the k-mer, but did
     * happen when {@code refseq.updateScope=otherTaxaOnly} withheld the project's own genomes from
     * that update. OTHER is the right answer there too - the k-mer is carried by something that is
     * not one of these children - and it keeps such a database's index bounded rather than letting
     * it grow by a factor of however many genomes went unseen.
     *
     * @param leaf       the leaf resolved for the current region, or {@code null} if there is none
     * @param storedNode the node the k-mer is stored at, i.e. the one being refined
     * @return the index the pair is to be recorded under
     */
    static int childIndexUnder(SmallTaxTree.SmallTaxIdNode leaf, SmallTaxTree.SmallTaxIdNode storedNode) {
        for (SmallTaxTree.SmallTaxIdNode node = leaf; node != null && node != storedNode; node = node.getParent()) {
            if (node.getParent() == storedNode) {
                return node.storeIndex;
            }
        }
        return OTHER_VALUE;
    }

    /**
     * Decides whether a k-mer read from a reference sequence takes part in this pass at all. Asked once
     * per k-mer, before it is looked up in the store, so that a goal which does not need every k-mer
     * pays neither for the lookup nor for what follows it.
     * <p>
     * Every k-mer takes part unless a subclass says otherwise. One that answers selectively has to do
     * so as a function of the k-mer alone, so that a k-mer met again - which is the normal case here,
     * as every genome of a taxon is read - is treated the same way every time.
     *
     * @param kmer the k-mer, encoded as a {@code long}
     * @return whether it should be looked up and recorded
     */
    protected boolean considerKMer(long kmer) {
        return true;
    }

    /**
     * Fills {@link #refinementBits} with the store index of every node of {@code statNodes} that
     * matches a configured {@link FTConfigKey.RefinementPosition}.
     *
     * @param statNodes the stored nodes to test, i.e. those carrying at least one k-mer
     * @throws IllegalStateException if the tax tree's store indices cannot be used as bit indices
     */
    private void initRefinementBits(Set<SmallTaxTree.SmallTaxIdNode> statNodes) {
        refinementBits = new BitSet();
        for (SmallTaxTree.SmallTaxIdNode s : statNodes) {
            if (s != null && FTConfigKey.RefinementPosition.getMatchingNodeFor(s, refinementIntervals) != null) {
                refinementBits.set(s.storeIndex);
            }
        }
    }

    /**
     * Returns whether the given stored node is selected for refinement, by testing its store index's
     * bit in {@link #refinementBits}. Indices beyond the highest marked one simply read as unset.
     *
     * @param node the stored node to test
     * @return whether the node is selected for refinement
     */
    protected boolean isRefinementNode(SmallTaxTree.SmallTaxIdNode node) {
        return refinementBits.get(node.storeIndex);
    }

    /**
     * Creates the {@link MyFastaReader} used to scan the reference fastas and record leaf-node
     * indices, configured from the current project's k-mer and genome-selection
     * settings.
     *
     * @param regionsPerTaxid trie of the fasta regions to read per taxon
     * @return the fasta reader to use for this pass
     */
    @Override
    protected AbstractStoreFastaReader createFastaReader(AbstractRefSeqFastaReader.StringLong2DigitTrie regionsPerTaxid) {
        MyFastaReader reader = new MyFastaReader(intConfigValue(GSConfigKey.FASTA_LINE_SIZE_BYTES),
                relevantNodes,
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
        readers.add(reader);
        return reader;
    }

    /**
     * Fasta reader that hands the enclosing goal, for each k-mer already contained in the k-mer store
     * below a taxon selected for refinement, the store index of the leaf tax node it stems from. It
     * relaxes the usual per-taxon limits so that the subsequent clustering phase can work with the
     * maximum amount of data.
     */
    protected class MyFastaReader extends AbstractUpdateFastaReader
            implements RadixKMerStore.BatchValueConsumer<SmallTaxTree.SmallTaxIdNode> {
        private long entries;
        /**
         * Buffers for the batched store lookup, or {@code null} if the k-mer store does not support
         * batching, in which case {@link #handleStore(long)} looks every k-mer up on its own. One instance
         * per reader, i.e. per thread, as the buffers are not thread-safe.
         */
        private final RadixKMerStore.BatchBuffers batch;

        /**
         * Creates the reader.
         *
         * @param bufferSize             size of the fasta read buffer in bytes
         * @param taxNodes               the tax nodes to be considered
         * @param accessionMap           accession-to-tax-node mapping, or {@code null} if unused
         * @param k                      the k-mer size
         * @param maxGenomesPerTaxId     maximum number of genomes to read per taxon
         * @param maxGenomesPerTaxIdRank rank at which the per-taxon genome limit applies
         * @param maxKmersPerTaxId       maximum number of k-mers to read per taxon
         * @param maxDust                maximum allowed low-complexity (dust) content
         * @param kMerSampling               k-mer sampling step size
         * @param assemblyAccessionsOnly    whether only complete genomes are considered
         * @param regionsPerTaxid        trie of the fasta regions to read per taxon
         * @param enableLowerCaseBases   whether lower-case bases are treated as valid
         */
        public MyFastaReader(int bufferSize, Set<TaxTree.TaxIdNode> taxNodes, AccessionMap accessionMap,
                             int k, int maxGenomesPerTaxId, Rank maxGenomesPerTaxIdRank, long maxKmersPerTaxId, int maxDust, int kMerSampling, boolean assemblyAccessionsOnly, StringLong2DigitTrie regionsPerTaxid, boolean enableLowerCaseBases) {
            super(bufferSize, taxNodes, accessionMap, k, maxGenomesPerTaxId, maxGenomesPerTaxIdRank, maxKmersPerTaxId, maxDust, kMerSampling, assemblyAccessionsOnly, regionsPerTaxid, enableLowerCaseBases, booleanConfigValue(GSConfigKey.ID_NODES), booleanConfigValue(GSConfigKey.FILE_NODES), booleanConfigValue(GSConfigKey.DATA_NODES));
            batch = kmerStore instanceof RadixKMerStore ? new RadixKMerStore.BatchBuffers(BATCH_SIZE) : null;
        }

        /**
         * Returns the compact tax tree of the database being indexed.
         *
         * @return the database's {@link SmallTaxTree}
         */
        @Override
        protected SmallTaxTree getTree() {
            return smallTaxTree;
        }

        // The following 3 override ensure that all k-mers for considers tax nodes are used
        // so that the clustering phase can work with the maximum amount of data.
        /**
         * Processes a fasta header line, resolving the current tax node and marking its region for
         * inclusion whenever the node is among the relevant nodes, so that all of its k-mers are
         * considered. Also refreshes the current leaf node.
         */
        @Override
        protected void infoLine() {
            if (ignoreMap) {
                node = mappedNode;
            } else {
                updateNodeFromInfoLine();
            }
            if (node != null && taxNodes.contains(node)) {
                includeRegion = true;
            }
            updateLeafNode();
        }

        /**
         * Does nothing; region ends are intentionally ignored so that no k-mers are dropped when a
         * region completes.
         */
        @Override
        protected void endRegion() {
            // Intentionally empty.
        }

        /**
         * Always allows further k-mers to be read, disabling the usual per-taxon cap so that the
         * clustering phase has the maximum amount of data available.
         *
         * @return always {@code true}
         */
        @Override
        public boolean isAllowMoreKmers() {
            return true;
        }

        /**
         * Registers the current k-mer in the filter. If the k-mer is stored in the database below a
         * taxon selected for refinement, the pair is recorded under the direct child of that taxon
         * the region's leaf lies under, or under {@link #OTHER_VALUE} when there is none - see
         * {@link #childIndexUnder}.
         *
         * @return {@code true} if a new (k-mer, index) entry was added, {@code false} otherwise
         */
        @Override
        protected boolean handleStore(long kmer) {
            if (!considerKMer(kmer)) {
                // Asked before the store is consulted, so that a goal which only estimates can skip the
                // lookup - which is the expensive part of this pass - and not just its own bookkeeping.
                return false;
            }
            int index = leafNode == null ? OTHER_VALUE : leafNode.storeIndex;
            if (batch != null) {
                // The leaf node is only known while reading, so it travels with the k-mer as the batch's
                // payload and the filter is written in accept() once the whole batch has been looked up.
                if (batch.add(kmer, index)) {
                    flushBatch();
                }
                // Whether the k-mer ends up in the filter is not known yet. The return value only feeds
                // kmersInRegion, which this reader does not use: it caps nothing (isAllowMoreKmers() is
                // always true) and endRegion(), its only other consumer, is overridden to do nothing.
                return false;
            }
            SmallTaxTree.SmallTaxIdNode storedNode = kmerStore.getLong(kmer, null);
            // Checking storedNode for not null is very important because otherwise to many
            // k-mers get entered - many more as estimated for filter size from above.
            if (storedNode != null && isRefinementNode(storedNode)) {
                // The leaf is at hand here, so the walk starts from the node itself and no index has
                // to be resolved back into one, as it does on the batched path.
                return putIndex(kmer, childIndexUnder(leafNode, storedNode));
            }
            return false;
        }


        /**
         * Looks up all buffered k-mers in one batch and records the ones stored below a taxon selected
         * for refinement. Called when the buffer is full and, for the trailing k-mers, once all fastas
         * have been read.
         */
        protected void flushBatch() {
            if (batch != null && !batch.isEmpty()) {
                ((RadixKMerStore<SmallTaxTree.SmallTaxIdNode>) kmerStore).getBatch(batch, this);
            }
        }

        /**
         * Records one k-mer of a flushed batch, which by construction is stored in the database. Only
         * called for k-mers that are present, so the {@code null} check of the unbatched path is
         * implicit here.
         *
         * @param kmer       the k-mer that was looked up
         * @param index      the store index of the leaf node the k-mer stems from, as buffered with it
         * @param storedNode the node the k-mer is stored against
         */
        @Override
        public void accept(long kmer, int index, SmallTaxTree.SmallTaxIdNode storedNode) {
            if (isRefinementNode(storedNode)) {
                putIndex(kmer, childIndexUnder(leafForIndex(index), storedNode));
            }
        }

        /**
         * Returns the leaf node a buffered k-mer was read under, from the store index that travelled
         * with it in the batch.
         * <p>
         * The batch outlives the region it was filled in - {@link #endRegion()} is a no-op here, so a
         * partly filled buffer carries into the next one - which is why the leaf cannot simply be read
         * off {@link #leafNode} at this point and has to be recovered from what was buffered. The
         * store's index map is a plain array, so this is one indexed read.
         *
         * A negative index answers {@code null} as well. {@code Database.initStoreIndices()} gives
         * every node of the tree an index, so a leaf without one would mean the tree and the store
         * had come apart - but the index map is a bare array, and reading it at -1 would end the
         * pass with an out-of-bounds throw rather than with the OTHER slot such a pair belongs in.
         *
         * @param index the store index buffered with the k-mer, or {@link #OTHER_VALUE} if its region
         *              resolved to no leaf at all
         * @return the leaf node, or {@code null} if there was none
         */
        private SmallTaxTree.SmallTaxIdNode leafForIndex(int index) {
            return index == OTHER_VALUE || index < 0 ? null : kmerStore.getValueForIndex(index);
        }

        /**
         * Hands the (k-mer, leaf index) pair to whatever the concrete goal does with it.
         *
         * @param kmer  the k-mer to record
         * @param index the store index of the leaf node the k-mer stems from
         * @return whether the pair was newly added
         */
        private boolean putIndex(long kmer, int index) {
            // What becomes of the pair is the concrete goal's business; both of them do it without
            // a lock, one by setting bits atomically, the other into a sketch of its own.
            long h = KMerIndexFilterHelper.combine(kmer, index);
            if (record(h)) {
                entries++;
                return true;
            }
            return false;
        }

        /**
         * Returns the number of entries this reader added to the filter.
         *
         * @return the number of (k-mer, child index) pairs this reader newly inserted
         */
        public long getEntries() {
            return entries;
        }
    }
}
