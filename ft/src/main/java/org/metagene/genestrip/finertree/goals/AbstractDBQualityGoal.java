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
import org.metagene.genestrip.finertree.FTGoalKey;
import org.metagene.genestrip.finertree.FTProject;
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
import org.metagene.genestrip.store.RadixKMerStore;
import org.metagene.genestrip.tax.Rank;
import org.metagene.genestrip.tax.SmallTaxTree;
import org.metagene.genestrip.tax.TaxTree;
import org.metagene.genestrip.tax.TaxNodeSelection;

import java.io.File;
import java.util.*;

/**
 * What the two database-quality goals have in common: reading the genomic fasta files back, resolving
 * each record to the leaf tax node the fill filed it at, and looking its *k*-mers up in the database.
 * Every *k*-mer the database also holds forms a (*k*-mer, leaf) pair.
 * <p>
 * What becomes of a pair is left to {@link MyFastaReader#count}, and that is where the two goals part:
 * {@link DBQualitySizeGoal} sketches the pairs to find out how many distinct ones there are, and
 * {@link DBQualityCountsGoal} deduplicates them through a filter of that size and tallies them.
 *
 * @param <T> the type of result the concrete goal produces
 * @param <P> the concrete FT project type
 */
public abstract class AbstractDBQualityGoal<T, P extends FTProject> extends FastaReaderGoal<T, P> implements Goal.LogHeapInfo {
    /**
     * Number of k-mers a reader buffers before looking them up in one batch. Sized to keep enough
     * independent cache-missing lookups in flight for the memory-level parallelism to pay, without the
     * per-batch bookkeeping outweighing it - the same value {@code AbstractKMerIndexGoal} uses.
     */
    protected static final int BATCH_SIZE = 128;

    /** Supplies the accession-to-tax-node map, or is left unasked where no RefSeq release is read. */
    protected final ObjectGoal<AccessionMap, P> accessionMapGoal;
    /** Supplies the database whose k-mers the genomes are compared against. */
    protected final ObjectGoal<Database, P> storeGoal;

    /** The database's taxonomy, held between {@link #prepare()} and {@link #releaseAfterPass()}. */
    protected SmallTaxTree tree;
    /** The database's k-mer store, held for the same stretch. */
    protected KMerStore<SmallTaxTree.SmallTaxIdNode> kMerSortedArray;
    /** The readers of the pass under way, one per thread, collected as they are created. */
    protected List<MyFastaReader> readers;
    /**
     * The tree's nodes by dense position.
     * <p>
     * The reading path addresses a node by its position and never by its tax id: the tally map of
     * {@link DBQualityCountsGoal} is keyed by a String, and looking a tally up in it twice per
     * (k-mer, data taxon) pair hashed a string a few billion times per run. {@link SmallTaxTree}
     * numbers its nodes densely for exactly this (see {@link SmallTaxTree#getNodeCount()}), so an
     * array does the same job with an indexed read.
     */
    protected SmallTaxTree.SmallTaxIdNode[] nodeByPos;
    /** Whether the node at that position is a leaf in the sense of {@link SmallTaxTree.SmallTaxIdNode#isLeaf}. */
    protected boolean[] leafByPos;
    /** The number of node positions, i.e. the length of {@link #nodeByPos} and {@link #leafByPos}. */
    protected int nodeCount;

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
    protected AbstractDBQualityGoal(P project, FTGoalKey key, ExecutionContext bundle, ObjectGoal<Set<RefSeqCategory>, P> categoriesGoal,
                                    ObjectGoal<TaxNodeSelection, P> taxNodesGoal,
                                    RefSeqFnaFilesDownloadGoal fnaFilesGoal,
                                    ObjectGoal<Map<File, TaxTree.TaxIdNode>, P> additionalGoal,
                                    ObjectGoal<AccessionMap, P> accessionMapGoal, ObjectGoal<Database, P> storeGoal,
                                    Goal<P>... deps) {
        super(project, key, bundle, categoriesGoal, taxNodesGoal, fnaFilesGoal, additionalGoal, Goal.append(deps, accessionMapGoal, storeGoal));
        this.storeGoal = storeGoal;
        this.accessionMapGoal = accessionMapGoal;
    }

    /**
     * Loads the database and works out, for every node of its tree, its dense position and whether it
     * is a leaf. Both goals need exactly this before they can read a single base.
     */
    protected void prepare() {
        GSProject project = getProject();
        // What a leaf is has to agree between the fill and this measure, and what secures that
        // agreement is an artificial node the fill gives to *every* taxon it files a genome at.
        // Then no taxonomy node holds a genome's k-mers, a node without children is exactly a node
        // the fill filed a genome at, and isLeaf() is one test that needs to know nothing about
        // ranks.
        //
        // `dataNodes' and `genomeNodes' both do that, and either will serve. ReworkingStoreFastaReader
        // .reworkNode() fires the DATA branch for any node not already of DATA rank, and the GENOME
        // branch for any node not already of GENOME rank, so with either on a taxon receiving a
        // contig always ends up with a child. `fileNodes' and `idNodes' do not qualify on their own:
        // a file node is created only where a file is in hand, and both are optional refinements
        // *below* whichever of the two above is in force.
        //
        // Requiring one of them is not what used to shut this goal out of a database refined below
        // the species -- that was the second half of the old guard, which REJECTED `fileNodes' while
        // only a DATA node could be recognised as a leaf. Both halves went at once; only the first
        // is coming back, now widened. `fileNodes' and `idNodes' stay free, and the `cdiff' project
        // of ft-db-exp2, which needs file nodes because the taxonomy supplies no children below the
        // species, has data nodes on as every project here does.
        if (!project.booleanConfigValue(GSConfigKey.DATA_NODES)
                && !project.booleanConfigValue(GSConfigKey.GENOME_NODES)) {
            throw new IllegalStateException(
                    "This goal requires data nodes (dataNodes=true) or genome nodes (genomeNodes=true)");
        }

        tree = storeGoal.get().getTaxTree();
        nodeCount = tree.getNodeCount();
        nodeByPos = new SmallTaxTree.SmallTaxIdNode[nodeCount];
        leafByPos = new boolean[nodeCount];
        for (SmallTaxTree.SmallTaxIdNode node : tree) {
            int pos = node.getPosition();
            nodeByPos[pos] = node;
            leafByPos[pos] = node.isLeaf();
        }
        kMerSortedArray = storeGoal.get().convertKMerStore();
        readers = new ArrayList<>();
    }

    /**
     * Releases what {@link #prepare()} loaded. Both goals call this from their own finally block, so
     * that an aborted pass leaves no database behind. What a pass allocated for itself -- a filter, a
     * sketch -- it releases itself, in the same block.
     */
    protected void releaseAfterPass() {
        tree = null;
        kMerSortedArray = null;
        readers = null;
        nodeByPos = null;
        leafByPos = null;
    }

    /**
     * Creates the fasta reader that compares genome *k*-mers against the database and registers it, so
     * that its last, partial batch can be flushed once the pass is done.
     *
     * @param contigsPerTaxid the trie counting contigs per tax id
     * @return the fasta reader to use for reading the genomic fasta files
     */
    @Override
    protected AbstractStoreFastaReader createFastaReader(AbstractRefSeqFastaReader.StringLong2DigitTrie contigsPerTaxid) {
        MyFastaReader reader = newReader(contigsPerTaxid);
        readers.add(reader);
        return reader;
    }

    /**
     * Creates the reader of this pass, which is the one thing the two goals do not share.
     *
     * @param contigsPerTaxid the trie counting contigs per tax id
     * @return the reader, configured from the project
     */
    protected abstract MyFastaReader newReader(AbstractRefSeqFastaReader.StringLong2DigitTrie contigsPerTaxid);

    /**
     * Fasta reader that, for each *k*-mer read from a genome, resolves the leaf tax node of the record
     * it was read in and checks whether the database holds the *k*-mer at all. Every *k*-mer it does
     * hold is handed to {@link #count}, which each goal answers for itself.
     */
    protected abstract class MyFastaReader extends AbstractUpdateFastaReader
            implements RadixKMerStore.BatchValueConsumer<SmallTaxTree.SmallTaxIdNode> {
        /** Number of (k-mer, leaf node) pairs this reader formed. */
        protected long entries;

        /**
         * Buffers for the batched store lookup, or {@code null} when it cannot be used. The store
         * lookup is memory-latency bound, and a batch lets many of its cache misses overlap - see
         * {@link RadixKMerStore#getBatch}.
         */
        private final RadixKMerStore.BatchBuffers batch;
        /** The current contig's leaf and its position, resolved once per contig rather than per k-mer. */
        private SmallTaxTree.SmallTaxIdNode cachedLeaf;
        private int cachedLeafPos = -1;

        /**
         * Creates the reader, reading everything but the contigs from the goal's configuration.
         *
         * @param contigsPerTaxid the trie counting contigs per tax id
         */
        protected MyFastaReader(StringLong2DigitTrie contigsPerTaxid) {
            super(intConfigValue(GSConfigKey.FASTA_LINE_SIZE_BYTES), taxNodesGoal.get().getSelected(),
                    isIncludeRefSeqFna() ? accessionMapGoal.get() : null,
                    intConfigValue(GSConfigKey.KMER_SIZE),
                    intConfigValue(GSConfigKey.MAX_DUST),
                    intConfigValue(GSConfigKey.KMER_SAMPLING),
                    booleanConfigValue(GSConfigKey.ASSEMBLY_ACCESSIONS_ONLY),
                    contigsPerTaxid,
                    booleanConfigValue(GSConfigKey.ENABLE_LOWERCASE_BASES),
                    booleanConfigValue(GSConfigKey.ID_NODES),
                    booleanConfigValue(GSConfigKey.GENOME_NODES),
                    booleanConfigValue(GSConfigKey.FILE_NODES),
                    booleanConfigValue(GSConfigKey.DATA_NODES));
            // Batched only while no per-taxon limit binds. A batched k-mer is counted after
            // handleStore() has already returned, so the return value can no longer say whether it
            // was, and that value feeds kmersInContig - which endContig() adds to the per-taxon
            // counters that maxGenomesPerTaxid is enforced from. At the default it does not bind and
            // nothing reads those counters back; with it set, the one-at-a-time path keeps the
            // accounting exact.
            boolean unlimited = intConfigValue(GSConfigKey.MAX_GENOMES_PER_TAXID) == Integer.MAX_VALUE;
            batch = (kMerSortedArray instanceof RadixKMerStore && unlimited)
                    ? new RadixKMerStore.BatchBuffers(BATCH_SIZE) : null;
        }

        /**
         * Looks the buffered k-mers up in one batch and hands those the database holds to
         * {@link #count}. Called when the buffer fills and, for the trailing ones, after all fastas
         * have been read.
         */
        protected void flushBatch() {
            if (batch != null && !batch.isEmpty()) {
                ((RadixKMerStore<SmallTaxTree.SmallTaxIdNode>) kMerSortedArray).getBatch(batch, this);
            }
        }

        /**
         * Resolves the contig's leaf as usual and caches its position, so that the reading path costs
         * a field read per k-mer instead of a lookup. The leaf changes per contig; the k-mers of a
         * contig are legion.
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
                        // worth a word, and once per contig rather than once per k-mer.
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
         * Looks the current *k*-mer up in the database and, if it holds it, forms the pair with the
         * leaf of the contig being read and hands it to {@link #count}. A *k*-mer read in a contig
         * that resolved to no leaf, or one the database does not hold, forms no pair.
         *
         * @param kmer the *k*-mer just read
         * @return what {@link #count} answered, or {@code false} where no pair was formed -- including
         * the batched case, where the pair is only formed once the batch comes back
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
         * Takes one k-mer of a flushed batch, which by construction the database holds - so the
         * {@code null} check of the unbatched path is implicit here - and hands it to {@link #count}.
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
         * Takes one (k-mer, leaf) pair, the k-mer being one the database holds.
         *
         * @param kmer       the k-mer
         * @param leafPos    the position of the leaf it was read in
         * @param storedNode the node the database stores it at
         * @return whether the pair was new, i.e. whether it was counted
         */
        protected abstract boolean count(long kmer, int leafPos, SmallTaxTree.SmallTaxIdNode storedNode);
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
}
