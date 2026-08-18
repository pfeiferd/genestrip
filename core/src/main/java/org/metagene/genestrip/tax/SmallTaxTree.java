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
package org.metagene.genestrip.tax;

import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;
import java.util.*;

import org.metagene.genestrip.tax.TaxTree.TaxIdNode;
import org.metagene.genestrip.util.ByteArrayUtil;
import org.metagene.genestrip.util.DigitTrie;

/**
 * A compact, serializable version of the taxonomy tree that only retains the nodes
 * required by a database and is used during read matching. In addition to tree
 * navigation and lowest-common-ancestor queries. Iterating the tree yields its nodes in
 * depth-first order.
 * <p>
 * The tree stays read-only while reads are matched: the per-read vote counters live in
 * thread-private arrays held by the matcher and indexed by the nodes' dense
 * {@link SmallTaxIdNode#getPosition()} (see {@link #reinitPositions()}). They used to sit in the
 * nodes with one slot per matcher thread, which put all threads' slots on one cache line and made
 * every vote on a hot node bounce that line between cores.
 */
public class SmallTaxTree implements Serializable, Iterable<SmallTaxTree.SmallTaxIdNode> {
	private static final long serialVersionUID = 1L;

	private transient Comparator<String> taxIdComparator;
	// Number of nodes, i.e. the number of positions reinitPositions() assigned. Transient like the
	// positions themselves, and set by every path that establishes them (build, load, restructure),
	// so that matching can read it without having to (re)assign anything.
	private transient int nodeCount;
	// Whether the node positions and depths currently describe this tree. Restructuring clears it and
	// reinitPositions() restores it; everything that depends on positions refuses to work in between,
	// so a tree that was changed cannot silently be matched against with stale numbering.
	private transient boolean positionsValid;
	private transient SmallTaxIdNode root;
	private transient DigitTrie<SmallTaxIdNode> taxIdNodeTrie;

	/**
	 * Builds a compact tree from the given full {@link TaxTree}, keeping only the nodes
	 * marked as required.
	 *
	 * @param taxTree the full taxonomy tree to derive the compact tree from
	 */
	public SmallTaxTree(TaxTree taxTree) {
		root = taxTree.getRoot() == null ? null : new SmallTaxIdNode(taxTree.getRoot());
		taxIdNodeTrie = new DigitTrie<SmallTaxIdNode>();
		root.initTrie(taxIdNodeTrie, 0);
		reinitPositions();
	}

	/**
	 * Serializes the tree by writing out its nodes starting from the root.
	 *
	 * @param out the stream to write the tree to
	 * @throws IOException if writing fails
	 */
	private void writeObject(ObjectOutputStream out) throws IOException {
		root.writeTree(out);
	}

	/**
	 * Deserializes the tree by reading its nodes and rebuilding the lookup trie.
	 *
	 * @param in the stream to read the tree from
	 * @throws IOException            if reading fails
	 * @throws ClassNotFoundException if a serialized class cannot be resolved
	 */
	private void readObject(ObjectInputStream in) throws IOException, ClassNotFoundException {
		root = SmallTaxIdNode.readTree(in);
		taxIdNodeTrie = new DigitTrie<SmallTaxIdNode>();
		root.initTrie(taxIdNodeTrie, 0);
		// Databases written before positions were densified carry the sparse numbering inherited from
		// the full taxonomy, so they are renumbered on load rather than having to be rebuilt. The
		// serialized format is unaffected - position is written either way.
		reinitPositions();
	}

	/**
	 * Replaces the sub-nodes of the node with the given tax id and re-registers the new
	 * descendants in the lookup trie.
	 *
	 * Restructuring invalidates the node positions, so {@link #reinitPositions()} has to run before
	 * the tree is used for matching again.
	 *
	 * @param taxId    the tax id of the node whose sub-nodes are replaced
	 * @param subNodes the new sub-nodes to set
	 * @return the modified node, or {@code null} if no node has the given tax id
	 */
	private SmallTaxIdNode setSubNodes(String taxId, SmallTaxIdNode[] subNodes) {
		SmallTaxIdNode node = taxIdNodeTrie.get(taxId);
		if (node != null) {
			node.setSubNodes(subNodes);
			// Re-register the replaced subtree and refresh its depths from this node's own depth.
			node.initTrie(taxIdNodeTrie, node.depth);
			// The tree gained or lost nodes, so the pre-order positions no longer describe it.
			positionsValid = false;
		}

		return node;
	}

	/**
	 * Replaces the sub-nodes of several nodes at once and re-establishes the node positions
	 * afterwards, so the tree is immediately usable again.
	 * <p>
	 * Restructuring is offered in bulk only, and renumbering is not the caller's job: the positions
	 * are the tree's own invariant, and a tree that has been changed but not renumbered must not
	 * escape. Renumbering once for the whole batch also keeps this linear in the number of nodes
	 * rather than linear per changed node.
	 *
	 * @param subNodesByTaxId the new sub-nodes per tax id; entries whose tax id is unknown are ignored
	 * @return the number of nodes whose sub-nodes were replaced
	 */
	public int setSubNodes(Map<String, SmallTaxIdNode[]> subNodesByTaxId) {
		int changed = 0;
		for (Map.Entry<String, SmallTaxIdNode[]> entry : subNodesByTaxId.entrySet()) {
			if (setSubNodes(entry.getKey(), entry.getValue()) != null) {
				changed++;
			}
		}
		reinitPositions();
		return changed;
	}

	/**
	 * Recomputes the {@code position} of every node by a fresh depth-first traversal, numbering them
	 * densely from zero, and returns the number of nodes.
	 * <p>
	 * Density is not cosmetic: matching keeps one counter per node and per consumer thread, so a
	 * sparse numbering would inflate those arrays - and their cache footprint - by whatever factor it
	 * is sparse. The nodes inherit their positions from the full taxonomy, of which this tree keeps
	 * only the required ones, so they arrive sparse (measured on the viral database: 46,000 nodes
	 * numbered up to 289,915). This is therefore called whenever a tree is built, loaded or
	 * restructured, and callers may rely on positions being dense. It is not public: the tree
	 * establishes its own numbering on construction, on load and after {@link #setSubNodes(Map)}, so
	 * no caller can be left holding a tree it has to renumber itself.
	 *
	 * @return the number of nodes in this tree
	 */
	private int reinitPositions() {
		nodeCount = root == null ? 0 : root.initPositions(0, 0) + 1;
		positionsValid = true;
		return nodeCount;
	}

	/**
	 * Returns whether the node positions and depths currently describe this tree, i.e. whether it has
	 * not been restructured since the last {@link #reinitPositions()}.
	 *
	 * @return whether the positions are valid
	 */
	public boolean isPositionsValid() {
		return positionsValid;
	}

	// Guards the one entry point every user of the positions has to pass first: a matcher sizes its
	// per-node arrays from getNodeCount() before it matches anything, so a restructured tree fails
	// loudly there instead of producing quietly wrong classifications. The per-node reads themselves
	// (getPosition(), and the depths isAncestorOf()/getLowestCommonAncestor() use) are deliberately
	// left unguarded: they run millions of times per file, while validity can only change between
	// runs, not during one.
	private void checkPositionsValid() {
		if (!positionsValid) {
			throw new IllegalStateException(
					"The tax tree was restructured without re-establishing its node positions.");
		}
	}

	/**
	 * Returns the number of nodes, which is also the number of distinct
	 * {@link SmallTaxIdNode#getPosition()} values - so an array indexed by position needs exactly this
	 * many entries.
	 * <p>
	 * It is established when the tree is built, loaded or restructured (see
	 * {@link #setSubNodes(Map)}), which is what lets matching size its per-node arrays without
	 * writing to the tree: the tree stays read-only while reads are matched.
	 *
	 * @return the number of nodes in this tree
	 */
	public int getNodeCount() {
		checkPositionsValid();
		return nodeCount;
	}


	/**
	 * Whether {@code ancestor} lies on the path from {@code node} up to the root,
	 * i.e. is an ancestor of {@code node} or {@code node} itself.
	 *
	 * @param node     the node whose ancestry is checked
	 * @param ancestor the candidate ancestor node
	 * @return {@code true} if {@code ancestor} is an ancestor of {@code node} or
	 *         {@code node} itself
	 */
	// Made final for potential inlining by JVM
	public final boolean isAncestorOf(SmallTaxIdNode node, final SmallTaxIdNode ancestor) {
		if (node == null || ancestor == null) {
			return false;
		}
		// The per-node depth (kept current by initTrie()/initPositions(), and relied on by
		// getLowestCommonAncestor() too) turns this into a bounded walk: an ancestor is never deeper
		// than its descendant, so a greater depth rules it out without touching memory at all, and
		// otherwise exactly depth-difference steps suffice. Walking to the root - the bulk of which
		// was wasted whenever the two nodes are unrelated - is what this replaces; that walk showed up
		// as the entire measurable cost of read classification in a profile of the matcher.
		int steps = node.depth - ancestor.depth;
		if (steps < 0) {
			return false;
		}
		while (steps > 0) {
			node = node.parent;
			steps--;
		}
		// == will do, faster than equals on works on closed set of nodes with equals not overriden.
		return node == ancestor;
	}

	/**
	 * Returns the lowest common ancestor of the two nodes, or {@code null} if they have
	 * none in common.
	 *
	 * @param node1 the first node
	 * @param node2 the second node
	 * @return the lowest common ancestor, or {@code null} if they have none in common
	 */
	// Made final for potential inlining by JVM
	public final SmallTaxIdNode getLowestCommonAncestor(final SmallTaxIdNode node1, final SmallTaxIdNode node2) {
		// Mild optimization
		if (node1 == node2) {
			return node1;
		}
		if (node1 == null || node2 == null) {
			return null;
		}
		// Align the deeper node to the shallower one, then walk both up in lock-step until they meet.
		// O(depth) using the per-node depth (kept current by initTrie()/initPositions()), and - unlike
		// scanning node2's whole ancestor chain for each of node1's ancestors, O(d1 * d2) - it stops as
		// soon as the paths join, so the common ancestor / nearby case (the bulk of look-ups) is cheap.
		SmallTaxIdNode a = node1;
		SmallTaxIdNode b = node2;
		while (a.depth > b.depth) {
			a = a.parent;
		}
		while (b.depth > a.depth) {
			b = b.parent;
		}
		while (a != b) {
			a = a.parent;
			b = b.parent;
		}
		// a == b now: their common ancestor, or null if the two nodes live in different trees.
		return a;
	}

	/**
	 * Sorts the given tax id strings in place by their nodes' position within the tree,
	 * falling back to lexicographic order for tax ids not present in the tree.
	 *
	 * @param taxids the list of tax id strings to sort
	 * @return the same list, sorted in place
	 */
	public List<String> sortTaxidsViaTree(List<String> taxids) {
		if (taxIdComparator == null) {
			taxIdComparator = new Comparator<String>() {
				@Override
				public int compare(String o1, String o2) {
					SmallTaxIdNode a = getNodeByTaxId(o1);
					SmallTaxIdNode b = getNodeByTaxId(o2);

					if (a == null && b == null) {
						return o1.compareTo(o2);
					}
					if (a == null) {
						return -1;
					}
					if (b == null) {
						return 1;
					}
					return a.compareTo(b);
				}
			};
		}
		Collections.sort(taxids, taxIdComparator);
		return taxids;
	}

	/**
	 * Returns the root node of the tree, or {@code null} if the tree is empty.
	 *
	 * @return the root node
	 */
	public SmallTaxIdNode getRoot() {
		return root;
	}

	/**
	 * Returns an iterator over all nodes of the tree in depth-first (pre-order) order.
	 */
	public Iterator<SmallTaxIdNode> iterator() {
		List<Integer> posL = new ArrayList<>();
		if (root != null) {
			posL.add(-1);
		}

		return new Iterator<SmallTaxIdNode>() {
			private SmallTaxIdNode nextNode = root;
			private List<Integer> posList = posL;

			@Override
			public SmallTaxIdNode next() {
				if (posL.size() == 0) {
					throw new NoSuchElementException();
				}
				SmallTaxIdNode res = nextNode;
				int nextPos = posList.get(posL.size() - 1) + 1;
				while (nextNode.subNodes == null || nextPos >= nextNode.subNodes.length) {
					posList.remove(posL.size() - 1);
					if (posL.isEmpty()) {
						break;
					}
					nextPos = posList.get(posL.size() - 1) + 1;
					nextNode = nextNode.parent;
				}
				if (!posL.isEmpty()) {
					nextNode = nextNode.subNodes[nextPos];
					posList.set(posL.size() - 1, nextPos);
					posList.add(-1);
				}
				return res;
			}

			@Override
			public boolean hasNext() {
				return posL.size() != 0;
			}
		};
	}


	/**
	 * Looks up the node with the given tax id via the lookup trie.
	 *
	 * @param taxId the tax id to look up
	 * @return the matching node, or {@code null} if none exists
	 */
	public SmallTaxIdNode getNodeByTaxId(String taxId) {
		return taxIdNodeTrie.get(taxId);
	}

	/**
	 * A node of a {@link SmallTaxTree}. Besides the usual tree links it holds the
	 * per-thread count slots used during matching and a {@code storeIndex} that links
	 * it to its entry in the k-mer database.
	 */
	public static class SmallTaxIdNode extends TaxIdInfo {
		private static final long serialVersionUID = 1L;

		/** The child nodes of this node, or {@code null} if it is a leaf. */
		private SmallTaxIdNode[] subNodes;
		/** Whether this node has been explicitly requested. */
		private boolean requested;

		/** The parent of this node, or {@code null} for the root. */
		protected transient SmallTaxIdNode parent;
		// Transient (not part of the serialized format, so old databases stay loadable) and reassigned
		// whenever the tree structure is established: initTrie() on build/load/attach and initPositions()
		// on reinit. Cached form of getLevel() (root = 0), used by getLowestCommonAncestor().
		/** The depth of this node (root = 0); see {@link #getLevel()}. */
		private transient int depth;
		// Made public for inlining
		/** Index linking this node to its entry in the k-mer database, or {@code -1} if unset. */
		public transient int storeIndex;

		/**
		 * Creates a node with the given tax id, name and rank and no sub-nodes.
		 *
		 * @param taxId the tax id of the node
		 * @param name  the name of the node
		 * @param rank  the taxonomic rank of the node
		 */
		public SmallTaxIdNode(String taxId, String name, Rank rank) {
			this(taxId, name, rank, null);
		}

		/**
		 * Creates a node with the given tax id, name, rank and sub-nodes, adopting the sub-nodes as its
		 * children.
		 * <p>
		 * Building a node complete is the supported way to give it children: the structure of a tree
		 * that is already in use may only be changed through {@link SmallTaxTree#setSubNodes(Map)},
		 * which re-establishes the node positions afterwards. This constructor is for nodes that are
		 * not attached to a tree yet.
		 *
		 * @param taxId    the tax id of the node
		 * @param name     the name of the node
		 * @param rank     the taxonomic rank of the node
		 * @param subNodes the sub-nodes to adopt, or {@code null} for a leaf
		 */
		public SmallTaxIdNode(String taxId, String name, Rank rank, SmallTaxIdNode[] subNodes) {
			super(taxId, rank);
			this.name = name;
			storeIndex = -1;
			if (subNodes == null) {
				this.subNodes = null;
			} else {
				setSubNodes(subNodes);
			}
		}

		private SmallTaxIdNode(TaxIdNode node) {
			super(node.taxId, node.rank);
			name = node.name;
			position = node.position;
			storeIndex = -1;
			List<TaxIdNode> tsubNodes = node.getSubNodes();
			if (tsubNodes != null) {
				int count = 0;
				for (int i = 0; i < tsubNodes.size(); i++) {
					TaxIdNode subNode = tsubNodes.get(i);
					if (subNode.isRequired()) {
						count++;
					}
				}
				if (count > 0) {
					subNodes = new SmallTaxIdNode[count];
					count = 0;
					for (int i = 0; i < tsubNodes.size(); i++) {
						TaxIdNode subNode = tsubNodes.get(i);
						if (subNode.isRequired()) {
							subNodes[count] = new SmallTaxIdNode(subNode);
							subNodes[count].parent = this;
							count++;
						}
					}
				}
			}
		}

		/**
		 * Returns the depth of this node in the tree (the root has level 0), computed by
		 * walking up to the root.
		 *
		 * @return the depth of this node in the tree
		 */
		// Only need in dbinfo - so, may be inefficient instead of field wasting some memory.
		public int getLevel() {
			int level = 0;
			for (SmallTaxIdNode current = parent; current != null; current = current.parent) {
				level++;
			}
			return level;
		}

		private int initPositions(int counter, int depth) {
			position = counter;
			this.depth = depth;
			if (subNodes != null) {
				for (int i = 0; i < subNodes.length; i++) {
					counter = subNodes[i].initPositions(counter + 1, depth + 1);
				}
			}
			return counter;
		}

		/**
		 * Returns whether this node has been explicitly requested.
		 *
		 * @return {@code true} if this node is requested
		 */
		public boolean isRequested() {
			return requested;
		}

		/**
		 * Sets whether this node has been explicitly requested.
		 *
		 * @param value the new requested flag
		 */
		public void setRequested(boolean value) {
			requested = value;
		}

		/**
		 * Returns the parent of this node, or {@code null} for the root.
		 *
		 * @return the parent node
		 */
		public SmallTaxIdNode getParent() {
			return parent;
		}

		/**
		 * Returns the index linking this node to its entry in the k-mer database.
		 *
		 * @return the store index, or {@code -1} if unset
		 */
		public final int getStoreIndex() {
			return storeIndex;
		}

		/**
		 * Sets the index linking this node to its entry in the k-mer database.
		 *
		 * @param storeIndex the new store index
		 */
		public void setStoreIndex(int storeIndex) {
			this.storeIndex = storeIndex;
		}

		/**
		 * Returns the child nodes of this node, or {@code null} if it is a leaf.
		 *
		 * @return the sub-nodes
		 */
		public SmallTaxIdNode[] getSubNodes() {
			return subNodes;
		}

		/**
		 * Returns the number of child nodes of this node.
		 *
		 * @return the number of sub-nodes, or {@code 0} if it is a leaf
		 */
		public int getNumberOfSubNodes() {
			return subNodes == null ? 0 : subNodes.length;
		}

		/**
		 * Sets this node's sub-nodes and updates each sub-node's parent link to this node.
		 * <p>
		 * Not public: restructuring a node that belongs to a tree would leave the tree's node positions
		 * describing a shape that no longer exists. Attached trees are restructured through
		 * {@link SmallTaxTree#setSubNodes(Map)}, which renumbers afterwards; a detached node is built
		 * complete via {@link #SmallTaxIdNode(String, String, Rank, SmallTaxIdNode[])}.
		 *
		 * @param subNodes the new sub-nodes to set
		 */
		void setSubNodes(SmallTaxIdNode[] subNodes) {
			this.subNodes = subNodes;
			for (SmallTaxIdNode subNode : subNodes) {
				subNode.parent = this;
			}
		}

		/**
		 * Returns the {@code DATA}-rank child of this node, also searching underneath any
		 * intervening {@code REFINED} nodes, or {@code null} if there is none.
		 *
		 * @return the {@code DATA}-rank descendant, or {@code null} if there is none
		 */
		public SmallTaxIdNode getDataChild() {
			if (subNodes == null) {
				return null;
			}
			for (int i = 0; i < subNodes.length; i++) {
				SmallTaxIdNode subNode = subNodes[i];
				if (subNode.getRankOrdinal() == Rank.DATA.ordinal()) {
					return subNode;
				}
				// There may have beenn refinements here as well,
				// so must search under any refined node but not below...
				if (subNode.getRank() == Rank.REFINED) {
					subNode = subNode.getDataChild();
					if (subNode != null) {
						return subNode;
					}
				}
			}
			return null;
		}

		/**
		 * Recursively searches this node's descendants for one whose name equals the given
		 * string, returning the first match or {@code null}.
		 *
		 * @param name the name to search for
		 * @return the first matching descendant, or {@code null} if none matches
		 */
		public SmallTaxIdNode getDescendantWithName(String name) {
			if (subNodes == null) {
				return null;
			}
			for (int i = 0; i < subNodes.length; i++) {
				SmallTaxIdNode descendant = subNodes[i];
				if (descendant.name.equals(name)) {
					return descendant;
				}
				descendant = descendant.getDescendantWithName(name);
				if (descendant != null) {
					return descendant;
				}
			}
			return null;
		}

		/**
		 * Recursively searches this node's descendants for one whose name equals the given
		 * byte sub-array range, returning the first match or {@code null}.
		 *
		 * @param array the byte array holding the name to search for
		 * @param start the start index (inclusive) of the name range
		 * @param end   the end index (exclusive) of the name range
		 * @return the first matching descendant, or {@code null} if none matches
		 */
		public SmallTaxIdNode getDescendantWithName(byte[] array, int start, int end) {
			if (subNodes == null) {
				return null;
			}
			for (int i = 0; i < subNodes.length; i++) {
				SmallTaxIdNode descendant = subNodes[i];
				if (ByteArrayUtil.equals(array, start, end, descendant.name)) {
					return descendant;
				}
				descendant = descendant.getDescendantWithName(array, start, end);
				if (descendant != null) {
					return descendant;
				}
			}
			return null;
		}

		private final void initTrie(DigitTrie<SmallTaxIdNode> trie, int depth) {
			this.depth = depth;
			trie.set(taxId, this);
			if (subNodes != null) {
				for (int i = 0; i < subNodes.length; i++) {
					subNodes[i].initTrie(trie, depth + 1);
				}
			}
		}

		private void writeTree(ObjectOutputStream out) throws IOException {
			out.writeShort(rank);
			out.writeUTF(taxId);
			out.writeUTF(name);
			out.writeInt(position);
			out.writeBoolean(requested);
			if (subNodes != null) {
				out.writeShort(subNodes.length);
				for (int i = 0; i < subNodes.length; i++) {
					subNodes[i].writeTree(out);
				}
			} else {
				out.writeShort(0);
			}
		}

		private static SmallTaxIdNode readTree(ObjectInputStream in) throws IOException, ClassNotFoundException {
			short rank = in.readShort();
			String taxId = in.readUTF();
			String name = in.readUTF();
			SmallTaxIdNode node = new SmallTaxIdNode(taxId, name, Rank.byOrdinal(rank));
			node.position = in.readInt();
			node.requested = in.readBoolean();
			int size = in.readShort();
			if (size > 0) {
				node.subNodes = new SmallTaxIdNode[size];
				for (int i = 0; i < size; i++) {
					node.subNodes[i] = readTree(in);
					node.subNodes[i].parent = node;
				}
			}
			return node;
		}
	}
}
