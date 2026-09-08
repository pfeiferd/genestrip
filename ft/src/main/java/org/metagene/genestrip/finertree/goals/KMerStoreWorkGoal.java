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

import me.tongfei.progressbar.ProgressBar;
import org.metagene.genestrip.GSConfigKey;
import org.metagene.genestrip.probfilter.ProbFilter;
import org.metagene.genestrip.finertree.FTConfigKey;
import org.metagene.genestrip.finertree.FTGoalKey;
import org.metagene.genestrip.finertree.FTProject;
import org.metagene.genestrip.finertree.probfilter.KMerIndexFilterHelper;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.store.Database;
import org.metagene.genestrip.store.KMerStore;
import org.metagene.genestrip.tax.SmallTaxTree;
import org.metagene.genestrip.tax.SmallTaxTree.SmallTaxIdNode;
import org.metagene.genestrip.util.progressbar.GSProgressBarCreator;
import org.metagene.genestrip.util.progressbar.GSProgressUpdate;

import java.util.List;

/**
 * Abstract base goal that visits every *k*-mer of a loaded database's *k*-mer store and, for each
 * taxonomy "parent" node that sits at a configured refinement rank, determines how the *k*-mer's
 * presence is spread across the parent's direct child subtrees. It loads the {@link Database} and the
 * associated {@link ProbFilter}, iterates the store, and for every *k*-mer computes a
 * {@code boolean[]} membership vector ({@code bits}) with one slot per child subtree plus a trailing
 * "OTHER" slot, testing membership via {@code checkChild} over the bloom filter. Concrete subclasses
 * plug in behaviour through the {@link #beforeKMerStoreWork()}, {@link #inKMerStoreWork} and
 * {@link #afterKMerStoreWork()} hooks.
 * <p>
 * That vector is the whole of what the filter is read for, here and in every subclass, which is why
 * the index pass records a pair under a direct child of the node being refined rather than under the
 * genome it came from.
 *
 * @param <T> the result type produced by the concrete goal
 * @param <P> the concrete FT project type
 */
public abstract class KMerStoreWorkGoal<T, P extends FTProject> extends ObjectGoal<T, P> {
    private static final int INITIAL_MAX_CHILDREN = 256;

    private final ObjectGoal<Database, P> storeGoal;
    /** The goal providing the k-mer index Bloom filter used for subtree-membership tests. */
    protected final ObjectGoal<ProbFilter, P> bloomFilterGoal;
    /** The refinement positions, i.e. the ranks and tax ids at which parent nodes are considered. */
    protected final List<FTConfigKey.RefinementPosition> refinementIntervals;
    /** The loaded database, available from {@link #beforeKMerStoreWork()} on. */
    protected Database database;
    /** The visited k-mer store of {@link #database}, available from {@link #beforeKMerStoreWork()} on. */
    protected KMerStore<SmallTaxIdNode> kMerStore;
    /** The loaded k-mer index Bloom filter, available from {@link #beforeKMerStoreWork()} on. */
    protected ProbFilter bloomFilter;
    /** Whether the node at a given {@link SmallTaxIdNode#getPosition() position} is refined. */
    private boolean[] refinementByPos;

    /**
     * Creates the goal, adding the store and bloom-filter goals as dependencies and reading the
     * configured refinement positions from the project configuration.
     *
     * @param project         the FT project
     * @param goalKey         the goal key identifying this goal
     * @param storeGoal       the goal providing the loaded database whose *k*-mer store is visited
     * @param bloomFilterGoal the goal providing the *k*-mer index bloom filter used for
     *                        subtree-membership tests
     * @param deps            further goals this goal depends on
     */
    @SafeVarargs
    public KMerStoreWorkGoal(P project, FTGoalKey goalKey, ObjectGoal<Database, P> storeGoal,
                             ObjectGoal<ProbFilter, P> bloomFilterGoal,
                             Goal<P>... deps) {
        super(project, goalKey, Goal.append(deps, storeGoal, bloomFilterGoal));
        this.storeGoal = storeGoal;
        this.bloomFilterGoal = bloomFilterGoal;
        refinementIntervals = (List<FTConfigKey.RefinementPosition>) configValue(FTConfigKey.REFINEMENT_POSITIONS);
    }

    /**
     * Cleans the underlying store goal, releasing the loaded database so that its memory can be
     * reclaimed once this goal no longer needs it.
     */
    protected void cleanStoreGoal() {
        storeGoal.cleanThis();
    }

    /**
     * Loads the database, its converted *k*-mer store and the bloom filter, then visits every stored
     * *k*-mer. For each *k*-mer whose indexed value is a parent node at a refinement rank, the child
     * membership vector {@code bits} (one slot per child subtree plus a trailing "OTHER" slot) and the
     * number of set slots ({@code spread}) are computed and passed to {@link #inKMerStoreWork}. The
     * {@link #beforeKMerStoreWork()} and {@link #afterKMerStoreWork()} hooks bracket the traversal, and
     * progress is reported via a progress bar.
     */
    @Override
    protected void doMakeThis() {
        try {
            database = storeGoal.get();
            kMerStore = database.convertKMerStore();
            bloomFilter = bloomFilterGoal.get();
            initRefinementByPos();

            beforeKMerStoreWork();

            long max = kMerStore.getEntries();
            long[] current = new long[1];
            GSProgressUpdate update = new GSProgressUpdate() {
                @Override
                public long current() {
                    return current[0];
                }

                @Override
                public long max() {
                    return max;
                }
            };
            try (ProgressBar pb = createProgressBar(update)) {
                kMerStore.visit(new KMerStore.IndexedKMerStoreVisitor<SmallTaxIdNode>() {
                    private boolean[] bits = new boolean[INITIAL_MAX_CHILDREN];
                    // The set slots of `bits', ascending, valid up to `spread'. Filled in the same
                    // sweep that fills `bits', so it costs nothing: the sweep is already as long as
                    // the child list and a consumer that wants only the set slots would otherwise
                    // have to scan the whole array again to find them.
                    private int[] setSlots = new int[INITIAL_MAX_CHILDREN];

                    @Override
                    public void nextValue(KMerStore<SmallTaxIdNode> store, long kmer, int index, long pos) {
                        current[0] = pos;
                        SmallTaxIdNode parent = store.getValueForIndex(index);
                        if (parent != null) {
                            SmallTaxIdNode[] children = getSubNodes(parent);
                            if (children != null && children.length > 0) {
                                int n;
                                for (n = bits.length; n <= children.length; n *= 2) {
                                }
                                if (n > bits.length) {
                                    bits = new boolean[n];
                                    setSlots = new int[n];
                                }
                                int spread = 0;
                                for (int i = 0; i < children.length; i++) {
                                    bits[i] = checkChild(children[i], kmer);
                                    if (bits[i]) {
                                        setSlots[spread++] = i;
                                    }
                                }
                                bits[children.length] = bloomFilter.containsLong(KMerIndexFilterHelper.combine(kmer, KMerIndexBloomGoal.OTHER_VALUE));
                                if (bits[children.length]) {
                                    setSlots[spread++] = children.length;
                                }
                                inKMerStoreWork(parent, pos, bits, spread, setSlots);
                            }
                        }
                    }

                    /**
                     * Returns whether any genome below the given child subtree carries the k-mer.
                     * <p>
                     * One probe, not a walk of the subtree: the index pass records a pair under the
                     * direct child of the node being refined and never deeper (see
                     * {@code AbstractKMerIndexGoal#childIndexUnder}), so the child's own slot already
                     * answers for everything below it. A filter written by a build older than that
                     * registered the deepest leaf instead and has to be regenerated - which such a
                     * database needs in any case, its filter having been sized by a bound that did
                     * not hold.
                     *
                     * @param node the direct child subtree to test
                     * @param kmer the k-mer being visited
                     * @return whether the subtree carries it
                     */
                    protected boolean checkChild(SmallTaxIdNode node, long kmer) {
                        return bloomFilter.containsLong(KMerIndexFilterHelper.combine(kmer, node.storeIndex));
                    }
                });
            }

            afterKMerStoreWork();
        } finally {
            database = null;
            kMerStore = null;
            bloomFilter = null;
            refinementByPos = null;
        }
    }

    /**
     * Returns the direct child nodes of the given parent, but only if the parent sits at one of the
     * configured refinement positions; otherwise returns {@code null} so that the parent is skipped.
     *
     * @param parent the candidate parent node
     * @return the parent's child nodes when it matches a refinement position, or {@code null} otherwise
     */
    private void initRefinementByPos() {
        SmallTaxTree tree = database.getTaxTree();
        refinementByPos = new boolean[tree.getNodeCount()];
        for (SmallTaxIdNode node : tree) {
            refinementByPos[node.getPosition()] =
                    FTConfigKey.RefinementPosition.getMatchingNodeFor(node, refinementIntervals) != null;
        }
    }

    /**
     * Returns the child subtrees whose membership is to be recorded for a k-mer stored at the given
     * node, or {@code null} where that node is not at a refinement position and its k-mers are of no
     * interest here.
     *
     * @param parent the node a k-mer is stored at
     * @return its child subtrees, or {@code null} if the node is not being refined
     */
    protected SmallTaxIdNode[] getSubNodes(SmallTaxIdNode parent) {
        // A bit test, and not the rank/tax id matching itself: whether a node sits at a refinement
        // position depends on the node alone, while this is asked once per stored k-mer - eighty
        // million times on a bacterial database. AbstractKMerIndexGoal precomputes the same thing for
        // the same reason.
        if (refinementByPos[parent.getPosition()]) {
            // Concrete children only: the OTHER child is the placeholder for what they do not
            // cover, so it gets the trailing slot and never one among them.
            return parent.getSubNodesWithoutOther();
        } else {
            return null;
        }
    }

    /**
     * Hook invoked once before the *k*-mer store traversal begins, allowing subclasses to initialise
     * any accumulators they need.
     */
    protected abstract void beforeKMerStoreWork();

    /**
     * Hook invoked for every visited *k*-mer whose indexed value is a parent node at a refinement rank.
     *
     * @param parent the parent taxonomy node the *k*-mer maps to
     * @param pos    the position of the *k*-mer within the store (also used for progress reporting)
     * @param bits   the child-membership vector: one slot per child subtree plus a trailing "OTHER"
     *               slot, each {@code true} when the *k*-mer occurs in that subtree
     * @param spread the number of {@code true} slots in {@code bits}, i.e. how many subtrees the
     *               *k*-mer is spread across
     * @param setSlots the indices of those {@code true} slots, filled up to {@code spread}, so that a
     *               subclass can walk them without scanning {@code bits}
     */
    protected abstract void inKMerStoreWork(SmallTaxIdNode parent, long pos, boolean[] bits, int spread,
                                            int[] setSlots);

    /**
     * Hook invoked once after the *k*-mer store traversal has completed, allowing subclasses to
     * finalise and publish their result.
     */
    protected abstract void afterKMerStoreWork();

    /**
     * Creates a progress bar tracking the *k*-mer store traversal, or {@code null} when progress-bar
     * output is disabled in the configuration.
     *
     * @param update the callback providing the current and maximum progress values
     * @return a new progress bar, or {@code null} if progress bars are disabled
     */
    protected ProgressBar createProgressBar(GSProgressUpdate update) {
        return booleanConfigValue(GSConfigKey.PROGRESS_BAR) ?
                GSProgressBarCreator.newGSProgressBar(getKey().getName(), update.max(), 1000, " kmers", update, null, true) :
                null;
    }
}
