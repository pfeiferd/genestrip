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

import org.metagene.genestrip.probfilter.ProbFilter;
import org.metagene.genestrip.finertree.FTGoalKey;
import org.metagene.genestrip.finertree.FTProject;
import org.metagene.genestrip.finertree.cluster.DendrogramNode;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.store.Database;
import org.metagene.genestrip.store.KMerStore;
import org.metagene.genestrip.tax.Rank;
import org.metagene.genestrip.tax.SmallTaxTree;
import org.metagene.genestrip.tax.TaxTree;

import java.util.*;

/**
 * Rewrites the database's k-mer store and its {@link SmallTaxTree} so that they reflect the refined
 * structure computed by the clustering phase. For every refined parent taxon the dendrogram is
 * turned into new {@link SmallTaxTree.SmallTaxIdNode}s with fresh store value indices, each k-mer is
 * reassigned to the most specific refined node it matches, and the tax tree's subnode structure and
 * positions are updated accordingly.
 *
 * @param <P> the concrete {@link FTProject} type this goal operates on
 */
public class UpdateStoreGoal<P extends FTProject> extends KMerStoreWorkGoal<Database, P> implements Goal.LogHeapInfo {
    private final ObjectGoal<Map<SmallTaxTree.SmallTaxIdNode, DendrogramNode>, P> dendrogramGoal;

    private int idCounter;
    private KMerStore<String> orgkMerSortedArray;
    private Map<SmallTaxTree.SmallTaxIdNode, DendrogramNode> dendrograms;
    private Map<String, BitSetsForNodes> parentToBitSets;

    /**
     * Creates the goal that updates the k-mer store and tax tree to reflect the refined dendrograms.
     *
     * @param project         the project this goal belongs to
     * @param storeGoal       goal supplying the database whose k-mer store and tax tree are updated
     * @param dendrogramGoal  goal supplying, per refined parent taxon, the root dendrogram node
     * @param bloomFilterGoal goal supplying the k-mer index bloom filter driving the per-k-mer work
     * @param deps            further goals this goal depends on
     */
    @SafeVarargs
    public UpdateStoreGoal(P project, ObjectGoal<Database, P> storeGoal, ObjectGoal<Map<SmallTaxTree.SmallTaxIdNode, DendrogramNode>, P> dendrogramGoal, ObjectGoal<ProbFilter, P> bloomFilterGoal, Goal<P>... deps) {
        super(project, FTGoalKey.UPDATE_STORE_GOAL, storeGoal, bloomFilterGoal, append(deps, dendrogramGoal));
        this.dendrogramGoal = dendrogramGoal;
    }

    @Override
    protected void doMakeThis() {
        try {
            super.doMakeThis();
        } finally {
            orgkMerSortedArray = null;
            dendrograms = null;
            parentToBitSets = null;
        }
    }

    /**
     * Prepares the update by capturing the database's current k-mer store and, for every refined
     * parent taxon whose dendrogram still needs expansion, building a {@link BitSetsForNodes} that
     * creates the new refined tax nodes and their k-mer matching bit sets.
     */
    @Override
    protected void beforeKMerStoreWork() {
        orgkMerSortedArray = database.getKmerStore();
        dendrograms = dendrogramGoal.get();
        parentToBitSets = new HashMap<>();
        for (SmallTaxTree.SmallTaxIdNode key : dendrograms.keySet()) {
            DendrogramNode root = dendrograms.get(key);
            if (root != null) {
                if (root.getValueIndex() == -1) {
                    parentToBitSets.put(key.getTaxId(), new BitSetsForNodes(key, root));
                } else {
                    // Nothing to do...
                }
            }
        }
    }

    /**
     * Reassigns a single k-mer to the most specific refined node it matches. Given the k-mer's
     * original parent taxon and the bit set of leaf nodes it occurs in, the best matching refined
     * node is looked up and, if found, the k-mer's value index in the store is set to that node's
     * store index.
     *
     * @param parent the original parent taxon of the k-mer
     * @param pos    the k-mer's position in the k-mer store
     * @param bits   the bit set of leaf nodes the k-mer occurs in
     * @param spread the number of leaf nodes the k-mer spreads across
     */
    @Override
    protected void inKMerStoreWork(SmallTaxTree.SmallTaxIdNode parent, long pos, boolean[] bits, int spread,
                                   int[] setSlots) {
        BitSetsForNodes bitSetsForNodes = parentToBitSets.get(parent.getTaxId());
        if (bitSetsForNodes != null) {
            SmallTaxTree.SmallTaxIdNode node = bitSetsForNodes.getBestMatchingNode(bits);
            if (node != null) {
                orgkMerSortedArray.setIndexAtPosition(pos, node.getStoreIndex());
            }
        }
    }

    /**
     * Finalizes the update after all k-mers have been reassigned: attaches the two new child nodes
     * of every refined parent to the tax tree, reinitializes the tree's node positions, publishes the
     * updated database as this goal's result, and invalidates the original store goal whose content is
     * now stale.
     */
    @Override
    protected void afterKMerStoreWork() {
        SmallTaxTree tree = database.getTaxTree();
        // Adjust the small tree at each parent node now:
        // (It must be done later, cause the original tree is still needed in inKMerStoreWork().)
        Map<String, SmallTaxTree.SmallTaxIdNode[]> newSubNodesByTaxId = new HashMap<>();
        for (String key : parentToBitSets.keySet()) {
            BitSetsForNodes bitSetsForNodes = parentToBitSets.get(key);
            newSubNodesByTaxId.put(key,
                    new SmallTaxTree.SmallTaxIdNode[] { bitSetsForNodes.child1, bitSetsForNodes.child2 });
        }
        // Applied in one go: the tree re-establishes its node positions itself once the batch is in.
        tree.setSubNodes(newSubNodesByTaxId);
        // The tree gained nodes, so the newly inserted ones still need their store value index.
        database.initStoreIndices();
        // inKMerStoreWork() reassigned k-mer values; the store's per-taxid counts are recomputed on read
        // and rebaked when the updated database is serialized (see AbstractKMerStore#writeObject), so the
        // stale counts from the loaded store are never observed - nothing to fix up here.

        set(database);

        // We have changed the original database, so that the original store goal's
        // content becomes invalid:
        cleanStoreGoal();
    }


    /**
     * Holds, for one refined parent taxon, the new refined tax nodes together with the bit sets that
     * describe which of the parent's original subnodes each refined node covers. The bit sets are
     * kept sorted by ascending cardinality so that {@link #getBestMatchingNode(boolean[])} returns
     * the most specific matching refined node.
     */
    /**
     * Whether every slot set in {@code contained} is also set in {@code container}, over the first
     * {@code containedLength} slots.
     * <p>
     * A refined node may take a k-mer only if it covers <em>every</em> child the k-mer was found in.
     * The last of those slots is the OTHER one, which says the k-mer was read from a genome that is
     * not below any child at all - the node being refined itself, or something outside it. No refined
     * node covers that, so a k-mer carrying it belongs where it is and must not be pushed down: doing
     * so would claim a specificity the genome in the OTHER slot contradicts, which is the very thing
     * that slot is computed for. Since the bit sets are only as long as the child list, the OTHER slot
     * lies past their end, and a comparison bounded by the container's length silently ignores it.
     * <p>
     * The bound is passed in rather than taken from {@code contained.length}: that array is the
     * visitor's buffer, sized to a power of two and reused for every parent, so its tail holds
     * whatever the last, possibly larger, parent left there.
     *
     * @param container       the slots a refined node covers, one per child of the node being refined
     * @param contained       the slots the k-mer was found in, one per child plus OTHER
     * @param containedLength how many slots of {@code contained} belong to this parent
     * @return whether the refined node covers all of them
     */
    static boolean coversAll(boolean[] container, boolean[] contained, int containedLength) {
        for (int i = 0; i < containedLength; i++) {
            if (contained[i] && (i >= container.length || !container[i])) {
                return false;
            }
        }
        return true;
    }

    private class BitSetsForNodes {
        private final SmallTaxTree.SmallTaxIdNode[] orgSubnodes;

        private final SmallTaxTree.SmallTaxIdNode parent;
        private final boolean[][] bitSets;
        private final SmallTaxTree.SmallTaxIdNode[] nodes;
        private int bitsetPosCounter;
        private final SmallTaxTree.SmallTaxIdNode child1;
        private final SmallTaxTree.SmallTaxIdNode child2;

        /**
         * Builds the refined tax nodes and matching bit sets for the given parent taxon from its
         * dendrogram. The two direct children of the dendrogram root become this parent's new
         * subnodes, all internal dendrogram nodes get fresh refined tax nodes and store value
         * indices, the corresponding bit sets are initialized, and finally the bit sets are sorted.
         *
         * @param parent the original parent taxon being refined
         * @param root   the root dendrogram node describing the refinement of the parent
         */
        public BitSetsForNodes(SmallTaxTree.SmallTaxIdNode parent, DendrogramNode root) {
            // Concrete children only. Where the build already made an OTHER node, it is the
            // bucket and not one of the children it stands apart from; counting it here would
            // shift every slot and give it a second slot of its own.
            this.parent = parent;
            this.orgSubnodes = parent.getSubNodesWithoutOther();
            this.bitSets = new boolean[root.size() - 1 - orgSubnodes.length][];
            this.nodes = new SmallTaxTree.SmallTaxIdNode[bitSets.length];

            for (int i = 0; i < bitSets.length; i++) {
                bitSets[i] = new boolean[orgSubnodes.length];
            }
            bitsetPosCounter = 0;
            child1 = createNode(root.getChild1());
            child2 = createNode(root.getChild2());
            bitsetPosCounter = 0;
            initBitSets(root.getChild1());
            initBitSets(root.getChild2());
            sort();
        }

        /**
         * Sorts the refined nodes and their bit sets in place by ascending cardinality (number of
         * covered original subnodes), so that more specific nodes are tested first when matching.
         */
        public void sort() {
            // Very basic max sort is sufficient -
            // unfortunately, standard library methods don't work for this case.
            for (int i = 0; i < bitSets.length; i++) {
                int minIndex = i;
                int minCard = cardinality(bitSets[i]);
                for (int j = i + 1; j < bitSets.length; j++) {
                    int c = cardinality(bitSets[j]);
                    if (c < minCard) {
                        minIndex = j;
                        minCard = c;
                    }
                }
                boolean[] h = bitSets[i];
                bitSets[i] = bitSets[minIndex];
                bitSets[minIndex] = h;
                SmallTaxTree.SmallTaxIdNode hn = nodes[i];
                nodes[i] = nodes[minIndex];
                nodes[minIndex] = hn;
            }
        }

        private int cardinality(boolean[] bits) {
            int cardinality = 0;
            for (int i = 0; i < bits.length; i++) {
                if (bits[i]) {
                    cardinality++;
                }
            }
            return cardinality;
        }

        /**
         * Returns the most specific refined node whose covered subnode set contains all leaf nodes
         * given by {@code bits}. Because the bit sets are sorted by ascending cardinality, the first
         * containing node encountered is the smallest (most specific) match.
         *
         * @param bits the bit set of leaf nodes a k-mer occurs in
         * @return the best matching refined node, or {@code null} if none contains the given bits
         */
        public SmallTaxTree.SmallTaxIdNode getBestMatchingNode(boolean[] bits) {
            // orgSubnodes.length + 1, because that is how many slots of `bits' belong to this parent:
            // one per child and the trailing OTHER. The array itself is the visitor's buffer, grown to
            // a power of two and reused across parents, so anything past that is another parent's
            // leftovers and must not be read.
            for (int i = 0; i < bitSets.length; i++) {
                if (coversAll(bitSets[i], bits, orgSubnodes.length + 1)) {
                    return nodes[i];
                }
            }
            return null;
        }



        /**
         * Recursively turns a dendrogram node into a tax node. Internal dendrogram nodes (and the
         * "OTHER" bucket) get a freshly generated tax id, a new store value index, a generated name,
         * the {@link Rank#REFINED} rank, and their two children created recursively; leaf dendrogram
         * nodes are mapped back to the corresponding original subnode of the parent.
         *
         * @param node the dendrogram node to convert
         * @return the created or resolved {@link SmallTaxTree.SmallTaxIdNode}
         */
        protected SmallTaxTree.SmallTaxIdNode createNode(DendrogramNode node) {
            int valueIndex = node.getValueIndex();
            if (valueIndex == orgSubnodes.length) {
                // The OTHER bucket. Where the build already made an OTHER node under this parent,
                // that node is this bucket: making another would leave two of them under one
                // unrefined parent, meaning the same thing. It still needs a store value index, which
                // it does not carry from the build - nothing was filed there.
                SmallTaxTree.SmallTaxIdNode existing = parent.getOtherChild();
                if (existing != null) {
                    existing.setStoreIndex(orgkMerSortedArray.getAddValueIndex(existing.getTaxId()));
                    return existing;
                }
            }
            if (valueIndex == -1 || valueIndex == orgSubnodes.length) {
                String taxId = "000" + idCounter++;
                int index = orgkMerSortedArray.getAddValueIndex(taxId);
                StringBuilder name = new StringBuilder();
                buildName(node, name);
                SmallTaxTree.SmallTaxIdNode newNode;
                if (node.getValueIndex() == -1) {
                    // The slot is taken before the children are built so that they keep the positions
                    // they had when the parent was registered first; only the node object itself is
                    // created later, because a node receives its children through its constructor.
                    int bitsetPos = bitsetPosCounter++;
                    SmallTaxTree.SmallTaxIdNode[] newSubnodes = new SmallTaxTree.SmallTaxIdNode[2];
                    newSubnodes[0] = createNode(node.getChild1());
                    newSubnodes[1] = createNode(node.getChild2());
                    newNode = new SmallTaxTree.SmallTaxIdNode(taxId, name.toString(), Rank.REFINED, newSubnodes);
                    nodes[bitsetPos] = newNode;
                } else {
                    newNode = new SmallTaxTree.SmallTaxIdNode(taxId, name.toString(), Rank.REFINED);
                }
                newNode.setStoreIndex(index);
                return newNode;
            } else {
                return orgSubnodes[valueIndex];
            }
        }

        /**
         * Builds a human-readable name for a refined node by combining the tax ids of the first and
         * last original subnodes it covers, separated by {@code "/"} for adjacent leaves or
         * {@code "/.../"} when further leaves lie in between.
         *
         * @param node the dendrogram node to name
         * @param name the buffer the generated name is appended to
         */
        protected void buildName(DendrogramNode node, StringBuilder name) {
            int depth = buildFirstLastName(node, true, name);
            if (depth == 1) {
                name.append('/');
            }
            else if (depth > 1) {
                name.append("/.../");
            }
            if (depth > 0) {
                buildFirstLastName(node, false, name);
            }
        }

        /**
         * Descends the dendrogram along its first (leftmost) or last (rightmost) children until a
         * leaf is reached and appends that leaf's tax id (or {@code "OTHER"}) to the buffer,
         * returning the number of internal nodes traversed on the way.
         *
         * @param node  the dendrogram node to descend from
         * @param first {@code true} to follow first children, {@code false} to follow last children
         * @param name  the buffer the resolved leaf tax id is appended to
         * @return the depth, i.e. the number of internal nodes traversed to reach the leaf
         */
        protected int buildFirstLastName(DendrogramNode node, boolean first, StringBuilder name) {
            int index = node.getValueIndex();
            if (index != -1) {
                if (index < orgSubnodes.length) {
                    name.append(orgSubnodes[index].getTaxId());
                }
                else {
                    // The same word core builds an OTHER node's name from, taken from the one
                    // place it is written down (see TaxTree.otherNodeName).
                    name.append(Rank.OTHER.getName());
                }
                return 0;
            }
            else {
                return buildFirstLastName(first ? node.getChild1() : node.getChild2(), first, name) + 1;
            }
        }

        /**
         * Recursively fills the bit set of each internal dendrogram node with the union of the
         * original subnodes covered by its two children. Leaf dendrogram nodes are encoded as the
         * negative-offset return value {@code -valueIndex - 1} rather than a bit set index.
         *
         * @param node the dendrogram node to initialize bit sets for
         * @return the bit set index of an internal node, or the encoded leaf offset for a leaf node
         */
        protected int initBitSets(DendrogramNode node) {
            int valueIndex = node.getValueIndex();
            if (node.getValueIndex() == -1) {
                int res = bitsetPosCounter++;
                int a = initBitSets(node.getChild1());
                int b = initBitSets(node.getChild2());
                boolean[] target = bitSets[res];
                for (int i = 0; i < target.length; i++) {
                    target[i] = ((a < 0) ? (i == -a - 1) : bitSets[a][i]) || ((b < 0) ? (i == -b - 1) : bitSets[b][i]);
                }
                return res;
            } else {
                return -valueIndex - 1;
            }
        }
    }
}
