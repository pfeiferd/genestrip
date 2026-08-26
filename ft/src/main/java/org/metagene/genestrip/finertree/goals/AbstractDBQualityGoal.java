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
    /** Whether the node at that position is a leaf in the sense of {@link #isLeafNode}. */
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

        tree = storeGoal.get().getTaxTree();
        nodeCount = tree.getNodeCount();
        nodeByPos = new SmallTaxTree.SmallTaxIdNode[nodeCount];
        leafByPos = new boolean[nodeCount];
        for (SmallTaxTree.SmallTaxIdNode node : tree) {
            int pos = node.getPosition();
            nodeByPos[pos] = node;
            leafByPos[pos] = isLeafNode(node);
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
     * @param admittedGenomes the shared set of genome keys admitted so far
     * @return the fasta reader to use for reading the genomic fasta files
     */
    @Override
    protected AbstractStoreFastaReader createFastaReader(AbstractRefSeqFastaReader.StringLong2DigitTrie contigsPerTaxid, AbstractRefSeqFastaReader.GenomeKeyTrie admittedGenomes) {
        MyFastaReader reader = newReader(contigsPerTaxid, admittedGenomes);
        readers.add(reader);
        return reader;
    }

    /**
     * Creates the reader of this pass, which is the one thing the two goals do not share.
     *
     * @param contigsPerTaxid the trie counting contigs per tax id
     * @param admittedGenomes the shared set of genome keys admitted so far
     * @return the reader, configured from the project
     */
    protected abstract MyFastaReader newReader(AbstractRefSeqFastaReader.StringLong2DigitTrie contigsPerTaxid, AbstractRefSeqFastaReader.GenomeKeyTrie admittedGenomes);

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
         * @param admittedGenomes the shared set of genome keys admitted so far
         */
        protected MyFastaReader(StringLong2DigitTrie contigsPerTaxid, AbstractRefSeqFastaReader.GenomeKeyTrie admittedGenomes) {
            super(intConfigValue(GSConfigKey.FASTA_LINE_SIZE_BYTES), taxNodesGoal.get().getSelected(),
                    isIncludeRefSeqFna() ? accessionMapGoal.get() : null,
                    intConfigValue(GSConfigKey.KMER_SIZE),
                    intConfigValue(GSConfigKey.MAX_GENOMES_PER_TAXID),
                    (Rank) configValue(GSConfigKey.MAX_PER_TAXID_RANK),
                    longConfigValue(GSConfigKey.MAX_KMERS_PER_TAXID),
                    intConfigValue(GSConfigKey.MAX_DUST),
                    intConfigValue(GSConfigKey.KMER_SAMPLING),
                    booleanConfigValue(GSConfigKey.ASSEMBLY_ACCESSIONS_ONLY),
                    contigsPerTaxid, admittedGenomes,
                    booleanConfigValue(GSConfigKey.ENABLE_LOWERCASE_BASES),
                    booleanConfigValue(GSConfigKey.ID_NODES),
                    booleanConfigValue(GSConfigKey.FILE_NODES),
                    booleanConfigValue(GSConfigKey.DATA_NODES));
            // Batched only while no per-taxon limit binds. A batched k-mer is counted after
            // handleStore() has already returned, so the return value can no longer say whether it
            // was, and that value feeds kmersInContig - which endContig() adds to the per-taxon
            // counters that maxGenomesPerTaxid and maxKMersPerTaxid are enforced from. At the defaults
            // neither binds and nothing reads those counters back; with either set, the one-at-a-time
            // path keeps the accounting exact.
            boolean unlimited = intConfigValue(GSConfigKey.MAX_GENOMES_PER_TAXID) == Integer.MAX_VALUE
                    && longConfigValue(GSConfigKey.MAX_KMERS_PER_TAXID) == Long.MAX_VALUE;
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
     * Whether the given node is where a genomic file's k-mers come to rest, and therefore the unit
     * the measures of these two goals are taken over.
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
     * {@code DBQualityCountsGoal.CountingReader.count} catches and turns into an
     * {@link IllegalStateException}.
     * <p>
     * {@link Rank#REFINED} is deliberately not an origin rank. A refined node is inserted by the
     * refinement <em>above</em> the origin nodes and holds the k-mers it moved down there, so it is
     * internal in exactly the way a taxonomy node is, and the measures restricted to what sits above
     * the data (see {@code DBQualityCountsGoal.Counts.aggregateSubtree}) have to keep counting it.
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
}
