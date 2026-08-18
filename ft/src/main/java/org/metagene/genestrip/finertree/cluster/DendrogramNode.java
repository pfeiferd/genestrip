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
package org.metagene.genestrip.finertree.cluster;

/**
 * A node of a binary dendrogram produced by hierarchical agglomerative clustering. A node is either
 * a leaf, which carries the index of a single clustered item ({@link #getValueIndex()}), or an
 * internal node, which joins two child nodes at a given merge {@link #getSimilarity() similarity}.
 * Each node may additionally carry an arbitrary application-specific {@link #getValue() value}.
 */
public class DendrogramNode {
    /**
     * Callback interface for a depth-first, pre-/post-order traversal of a dendrogram driven by
     * {@link DendrogramNode#visit(Visitor)}.
     */
    public interface Visitor {
        /**
         * Invoked when the traversal enters the given node, before its children are visited.
         *
         * @param node the node being entered
         */
        public void preNode(DendrogramNode node);

        /**
         * Invoked when the traversal leaves the given node, after its children have been visited.
         *
         * @param node the node being left
         */
        public void postNode(DendrogramNode node);
    }

    private final int valueIndex;
    private final DendrogramNode child1;
    private final DendrogramNode child2;
    private double similarity;
    private Object value;

    /**
     * Creates an internal node joining two child nodes at the given merge similarity. Its value
     * index is set to {@code -1} to mark it as a non-leaf node.
     *
     * @param child1     the first child node
     * @param child2     the second child node
     * @param similarity the similarity at which the two children were merged
     */
    public DendrogramNode(DendrogramNode child1, DendrogramNode child2, double similarity) {
        this.valueIndex = -1;
        this.child1 = child1;
        this.child2 = child2;
        this.similarity = similarity;
    }

    /**
     * Creates a leaf node for a single clustered item identified by its index.
     *
     * @param valueIndex the index of the clustered item represented by this leaf
     * @param similarity the similarity value associated with the leaf (typically the item's
     *                   self-similarity)
     */
    public DendrogramNode(int valueIndex, double similarity) {
        this.valueIndex = valueIndex;
        this.child1 = null;
        this.child2 = null;
        this.similarity = similarity;
    }

    /**
     * Returns the index of the clustered item represented by this leaf, or {@code -1} if this node
     * is an internal (non-leaf) node.
     *
     * @return the item index of a leaf, or {@code -1} for an internal node
     */
    public int getValueIndex() {
        return valueIndex;
    }

    /**
     * Returns the first child of this internal node, or {@code null} if this node is a leaf.
     *
     * @return the first child node, or {@code null} for a leaf
     */
    public DendrogramNode getChild1() {
        return child1;
    }

    /**
     * Returns the second child of this internal node, or {@code null} if this node is a leaf.
     *
     * @return the second child node, or {@code null} for a leaf
     */
    public DendrogramNode getChild2() {
        return child2;
    }

    /**
     * Returns the similarity associated with this node: for an internal node this is the similarity
     * at which its two children were merged.
     *
     * @return the similarity value of this node
     */
    public double getSimilarity() {
        return similarity;
    }

    /**
     * Attaches an arbitrary application-specific value to this node.
     *
     * @param value the value to associate with this node
     */
    public void setValue(Object value) {
        this.value = value;
    }

    /**
     * Returns the application-specific value previously attached to this node, or {@code null} if
     * none was set.
     *
     * @return the value associated with this node, or {@code null}
     */
    public Object getValue() {
        return value;
    }

    /**
     * Traverses this node and its subtree depth-first, invoking {@link Visitor#preNode(DendrogramNode)}
     * when entering each node and {@link Visitor#postNode(DendrogramNode)} when leaving it. Children
     * are only descended into for internal nodes.
     *
     * @param visitor the visitor to notify for each visited node
     */
    public void visit(Visitor visitor) {
        visitor.preNode(this);
        if (valueIndex == - 1) {
            child1.visit(visitor);
            child2.visit(visitor);
        }
        visitor.postNode(this);
    }

    /**
     * Returns the total number of nodes in the subtree rooted at this node, counting both internal
     * nodes and leaves.
     *
     * @return the number of nodes in this subtree (1 for a leaf)
     */
    public int size() {
        if (valueIndex == - 1) {
            return 1 + child1.size() + child2.size();
        }
        else {
            return 1;
        }
    }
}
