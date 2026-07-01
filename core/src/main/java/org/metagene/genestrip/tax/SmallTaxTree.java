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
 * navigation and lowest-common-ancestor queries, it provides per-node counters that
 * can be updated concurrently by several matcher threads (one counter slot per thread)
 * and summed along the path from a node to the root. Iterating the tree yields its
 * nodes in depth-first order.
 */
public class SmallTaxTree implements Serializable, Iterable<SmallTaxTree.SmallTaxIdNode> {
	private static final long serialVersionUID = 1L;

	private transient Comparator<String> taxIdComparator;
	private transient int countSize;
	private transient SmallTaxIdNode root;
	private transient DigitTrie<SmallTaxIdNode> taxIdNodeTrie;
	private transient Object owner;

	/**
	 * Builds a compact tree from the given full {@link TaxTree}, keeping only the nodes
	 * marked as required.
	 *
	 * @param taxTree the full taxonomy tree to derive the compact tree from
	 */
	public SmallTaxTree(TaxTree taxTree) {
		root = taxTree.getRoot() == null ? null : new SmallTaxIdNode(taxTree.getRoot());
		taxIdNodeTrie = new DigitTrie<SmallTaxIdNode>();
		root.initTrie(taxIdNodeTrie);
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
		root.initTrie(taxIdNodeTrie);
	}

	/**
	 * Replaces the sub-nodes of the node with the given tax id and re-registers the new
	 * descendants in the lookup trie.
	 *
	 * @param taxId    the tax id of the node whose sub-nodes are replaced
	 * @param subNodes the new sub-nodes to set
	 * @return the modified node, or {@code null} if no node has the given tax id
	 */
	public SmallTaxIdNode setSubNodes(String taxId, SmallTaxIdNode[] subNodes) {
		SmallTaxIdNode node = taxIdNodeTrie.get(taxId);
		if (node != null) {
			node.setSubNodes(subNodes);
			node.initTrie(taxIdNodeTrie);
		}

		return node;
	}

	/**
	 * Recomputes the {@code position} of every node by a fresh depth-first traversal.
	 */
	public void reinitPositions() {
		root.initPositions(0);
	}

	/**
	 * Sets the number of parallel counter slots per node (typically the number of
	 * matcher threads). Can only be set once.
	 *
	 * @param countSize the number of parallel counter slots per node
	 * @throws IllegalArgumentException if {@code countSize} is not positive
	 * @throws IllegalStateException    if it was already initialized to a different value
	 */
	public void initCountSize(int countSize) {
		if (countSize <= 0) {
			throw new IllegalArgumentException("Initialization size must be >= 0.");
		}
		if (this.countSize > 0 && countSize != this.countSize) {
			throw new IllegalStateException("Count count size can only be initialized once.");
		}
		this.countSize = countSize;
	}

	/**
	 * Resets all node counts and claims the tree for the given owner. The tree must be
	 * released by its current owner via {@link #releaseOwner()} before another owner can
	 * reset it.
	 *
	 * @param owner the object claiming ownership of the tree
	 * @throws IllegalArgumentException if the tree is still owned by a different owner
	 */
	public void resetCounts(Object owner) {
		if (owner == null) {
			throw new NullPointerException("owner must not be null");
		}
		if (this.owner == null) {
			this.owner = owner;
		} else if (this.owner != owner) {
			throw new IllegalArgumentException("Tax tree not released by previous owner: " + this.owner);
		}
		root.resetCounts();
	}

	/**
	 * Releases the current owner claimed via {@link #resetCounts(Object)}.
	 */
	public void releaseOwner() {
		this.owner = null;
	}

	/**
	 * Increments the counter of the given node in slot {@code index}. The {@code initKey}
	 * identifies the current read so that counters left over from a previous read are
	 * reset rather than accumulated.
	 *
	 * @param node    the node whose counter is incremented
	 * @param index   the counter slot to increment
	 * @param initKey the key identifying the current read
	 */
	// Made final for potential inlining by JVM
	public final void incCount(final SmallTaxIdNode node, final int index, final long initKey) {
		node.incCount(index, initKey, countSize);
	}

	/**
	 * Sums the counts in slot {@code index} (matching {@code initKey}) along the path
	 * from the given node up to the root.
	 *
	 * @param node    the node to start summing from
	 * @param index   the counter slot to sum
	 * @param initKey the key identifying the current read
	 * @return the sum of the matching counts from the node to the root
	 */
	// Made final for potential inlining by JVM
	public final int sumCounts(SmallTaxIdNode node, final int index, final long initKey) {
		int res = 0;
		while (node != null) {
			if (node.counts != null && node.countsInitKeys != null && node.countsInitKeys[index] == initKey) {
				res += node.counts[index];
			}
			node = node.parent;
		}
		return res;
	}

	/**
	 * Walks from the given node up to the root, accumulating the counts in slot
	 * {@code index} (matching {@code initKey}), and returns the first (lowest) node at
	 * which the running sum reaches {@code threshold}, or {@code null} if it never does.
	 *
	 * @param node      the node to start summing from
	 * @param index     the counter slot to sum
	 * @param initKey   the key identifying the current read
	 * @param threshold the running sum to reach
	 * @return the lowest node where the running sum reaches {@code threshold}, or
	 *         {@code null} if it never does
	 */
	// Made final for potential inlining by JVM
	public final SmallTaxIdNode lowestNodeWhereSumAboveThreshold(SmallTaxIdNode node, final int index, final long initKey, int threshold) {
		int res = 0;
		while (node != null) {
			if (node.counts != null && node.countsInitKeys != null && node.countsInitKeys[index] == initKey) {
				res += node.counts[index];
				if (res >= threshold) {
					return node;
				}
			}
			node = node.parent;
		}
		// TODO: Better return null here and handle null in calling code?
		return null;
	}

	/**
	 * Returns the number of parallel counter slots per node.
	 *
	 * @return the configured count size
	 */
	public int getCountSize() {
		return countSize;
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
		while (node != null) {
			// == will do, faster than equals on works on closed set of nodes with equals not overriden.
			if (node == ancestor) {
				return true;
			}
			node = node.parent;
		}

		return false;
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
		for (SmallTaxIdNode res = node1; res != null; res = res.parent) {
			for (SmallTaxIdNode ancestor2 = node2; ancestor2 != null; ancestor2 = ancestor2.parent) {
				if (res == ancestor2) {
					return res;
				}
			}
		}
		return null;
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

	/* Apparently not needed:
	protected int initPositions(int counter, SmallTaxIdNode taxIdNode) {
		taxIdNode.position = counter;
		SmallTaxIdNode[] subNodes = taxIdNode.subNodes;
		if (subNodes != null) {
			for (int i = 0; i < subNodes.length; i++) {
				counter = initPositions(counter + 1, subNodes[i]);
			}
		}
		return counter;
	}
	*/

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
		private transient int[] counts;
		private transient long[] countsInitKeys;
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
			super(taxId, rank);
			this.name = name;
			this.parent = parent;
			subNodes = null;
			storeIndex = -1;
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

		private int initPositions(int counter) {
			position = counter;
			if (subNodes != null) {
				for (int i = 0; i < subNodes.length; i++) {
					counter = subNodes[i].initPositions(counter + 1);
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
		 *
		 * @param subNodes the new sub-nodes to set
		 */
		public void setSubNodes(SmallTaxIdNode[] subNodes) {
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

		/**
		 * Increments this node's counter in slot {@code index}, lazily allocating the
		 * count arrays of the given {@code size}. The {@code initKey} identifies the
		 * current read so a stale count from a previous read is reset to one.
		 *
		 * @param index   the counter slot to increment
		 * @param initKey the key identifying the current read
		 * @param size    the size of the count arrays to allocate lazily
		 */
		public final void incCount(final int index, final long initKey, final int size) {
			if (counts == null || countsInitKeys == null) {
				synchronized (this) {
					if (counts == null) {
						counts = new int[size];
						countsInitKeys = new long[size];
					}
				}
			}
			if (countsInitKeys[index] == initKey) {
				counts[index]++;
			} else {
				countsInitKeys[index] = initKey;
				counts[index] = 1;
			}
		}

		private final void resetCounts() {
			if (countsInitKeys != null) {
				for (int i = 0; i < countsInitKeys.length; i++) {
					countsInitKeys[i] = -1;
				}
			}
			if (subNodes != null) {
				for (int i = 0; i < subNodes.length; i++) {
					subNodes[i].resetCounts();
				}
			}
		}

		private final void initTrie(DigitTrie<SmallTaxIdNode> trie) {
			trie.set(taxId, this);
			if (subNodes != null) {
				for (int i = 0; i < subNodes.length; i++) {
					subNodes[i].initTrie(trie);
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
