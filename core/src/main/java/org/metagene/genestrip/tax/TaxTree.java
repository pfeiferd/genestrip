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

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;

import me.tongfei.progressbar.ProgressBar;
import org.metagene.genestrip.io.BufferedLineReader;
import org.metagene.genestrip.io.StreamProvider;
import org.metagene.genestrip.io.StreamingFileResource;
import org.metagene.genestrip.io.StreamingResource;
import org.metagene.genestrip.util.ByteArrayUtil;
import org.metagene.genestrip.util.DigitTrie;
import org.metagene.genestrip.util.progressbar.GSProgressBarCreator;

/**
 * The NCBI taxonomy tree, built from the {@code nodes.dmp} and {@code names.dmp} dump
 * files. Nodes are looked up by tax id via an internal trie. The tree supports
 * ancestor and lowest-common-ancestor queries and can create artificial child nodes
 * (data, file and id nodes) used to track where k-mers originate from. It can be
 * converted into a compact {@link SmallTaxTree} for use during matching.
 */
public class TaxTree {
	private static final int MAX_LINE_SIZE = 4096;

	// Atomic so concurrent reader threads (each synchronizing on a different parent node when
	// creating data/file/id nodes) cannot lose updates and generate duplicate artificial ids.
	private final AtomicInteger nextArtCounter = new AtomicInteger(1);

	/** Comparator ordering nodes by tax id, with {@code null} sorting first. */
	public static final Comparator<TaxIdNode> NODE_COMPARATOR = new Comparator<TaxIdNode>() {
		@Override
		public int compare(TaxIdNode a, TaxIdNode b) {
			if (a == null && b == null) {
				return 0;
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

	/** File name of the NCBI taxonomy nodes dump. */
	public static final String NODES_DMP = "nodes.dmp";
	/** File name of the NCBI taxonomy names dump. */
	public static final String NAMES_DMP = "names.dmp";

	private final TaxIdNode root;
	private final TaxIdNodeTrie taxIdNodeTrie;

	/**
	 * Loads the taxonomy tree from the {@code nodes.dmp} and {@code names.dmp} files
	 * located in the given directory.
	 *
	 * @param path        directory containing the NCBI taxonomy dump files
	 * @param progressBar whether to show a progress bar while reading the files
	 */
	public TaxTree(File path, boolean progressBar) {
		try {
			taxIdNodeTrie = new TaxIdNodeTrie();
			File nodesFile = new File(path, NODES_DMP);
			File namesFile = new File(path, NAMES_DMP);
			if (progressBar) {
				try (StreamingResource.StreamAccess sa = new StreamingFileResource(nodesFile, true).openStream()) {
					try (ProgressBar pb = GSProgressBarCreator.newGSProgressBar("taxnodes", sa, null, true)) {
						root = readNodesFromStream(sa.getInputStream());
					}
				}
				try (StreamingResource.StreamAccess sa = new StreamingFileResource(namesFile, true).openStream()) {
					try (ProgressBar pb = GSProgressBarCreator.newGSProgressBar("taxnames", sa, null, true)) {
						readNamesFromStream(sa.getInputStream());
					}
				}
			}
			else {
				try (InputStream is = StreamProvider.getInputStreamForFile(nodesFile)) {
					root = readNodesFromStream(is);
				}
				try (InputStream is = StreamProvider.getInputStreamForFile(namesFile)) {
					readNamesFromStream(is);
				}
			}

			root.initPositions(0, 0);
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}
	
	/**
	 * Creates a compact, serializable {@link SmallTaxTree} from this tree.
	 *
	 * @return the newly created compact tax tree
	 */
	public SmallTaxTree toSmallTaxTree() {
		return new SmallTaxTree(this);
	}

	/**
	 * Whether {@code ancestor} lies on the path from {@code node} up to the root,
	 * i.e. is an ancestor of {@code node} or {@code node} itself.
	 *
	 * @param node     the node whose ancestry is checked
	 * @param ancestor the candidate ancestor node
	 * @return {@code true} if {@code ancestor} is {@code node} or one of its ancestors
	 */
	public boolean isAncestorOf(TaxIdNode node, TaxIdNode ancestor) {
		while (node != null) {
			if (node.equals(ancestor)) {
				return true;
			}
			node = node.parent;
		}

		return false;
	}

	/**
	 * Returns the lowest common ancestor of the two nodes, or {@code null} if they
	 * have none in common.
	 *
	 * @param node1 the first node
	 * @param node2 the second node
	 * @return the lowest common ancestor, or {@code null} if there is none
	 */
	public TaxIdNode getLowestCommonAncestor(final TaxIdNode node1, final TaxIdNode node2) {
		// Mild optimization
		if (node1 == node2) {
			return node1;
		}
		if (node1 == null || node2 == null) {
			return null;
		}
		// Align the deeper node to the shallower one, then walk both up in lock-step until they meet.
		// This is O(depth) using the per-node depth assigned by initPositions(), and - unlike scanning
		// node2's whole ancestor chain for each of node1's ancestors, O(d1 * d2) - it terminates as soon
		// as the paths join, so the common ancestor / nearby case (the bulk of update look-ups) is cheap.
		// Only parent pointers and the (immutable-after-build) depth are read, so it stays thread-safe.
		TaxIdNode a = node1;
		TaxIdNode b = node2;
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
	 * Reads a {@code names.dmp} stream and assigns names to the already loaded nodes,
	 * preferring scientific names.
	 *
	 * @param stream the {@code names.dmp} input stream to read
	 * @throws java.io.IOException if reading the stream fails
	 */
	protected void readNamesFromStream(InputStream stream) throws IOException {
		try (BufferedLineReader br = new BufferedLineReader(stream)) {
			int size;
			byte[] target = new byte[MAX_LINE_SIZE];

			while ((size = br.nextLine(target)) > 0) {
				boolean scientific = ByteArrayUtil.indexOf(target, 0, size, "scientific name") != -1;
				int a = ByteArrayUtil.indexOf(target, 0, size, '|');
				int b = ByteArrayUtil.indexOf(target, a + 1, size, '|');
				if (a != -1 && b != -1) {
					TaxIdNode node = taxIdNodeTrie.get(target, 0, a - 1);
					if (node != null) {
						String name = new String(target, a + 2, b - a - 3, StandardCharsets.UTF_8);
						if (node.name == null || scientific) {
							node.name = name;
						}
					}
				}
			}
		}
	}

	/**
	 * Reads a {@code nodes.dmp} stream, creating nodes, linking them to their parents
	 * and setting their ranks.
	 *
	 * @param stream the {@code nodes.dmp} input stream to read
	 * @return the root node (tax id "1")
	 * @throws java.io.IOException if reading the stream fails
	 */
	protected TaxIdNode readNodesFromStream(InputStream stream) throws IOException {
		TaxIdNode res = null;
		int size;
		byte[] target = new byte[MAX_LINE_SIZE];

		try (BufferedLineReader br = new BufferedLineReader(stream)) {
			while ((size = br.nextLine(target)) > 0) {
				int a = ByteArrayUtil.indexOf(target, 0, size, '|');
				int b = ByteArrayUtil.indexOf(target, a + 1, size, '|');
				int c = ByteArrayUtil.indexOf(target, b + 1, size, '|');
				if (a != -1 && b != -1) {
					TaxIdNode nodeA = taxIdNodeTrie.get(target, 0, a - 1, true);
					TaxIdNode nodeB = taxIdNodeTrie.get(target, a + 2, b - 1, true);

					if (nodeA != nodeB) {
						nodeB.addSubNode(nodeA);
					}
					Rank r = Rank.getRankFromBytes(target, b + 2, c - 1);
					nodeA.rank = (short) (r == null ? -1 : r.ordinal());

					if (nodeA == nodeB && "1".equals(nodeA.taxId)) {
						res = nodeA;
					}
				}
			}
		}

		return res;
	}

	/**
	 * Returns the artificial {@code DATA}-rank child of the given node, creating it
	 * (with a freshly generated artificial tax id) if it does not yet exist. Thread-safe.
	 *
	 * @param node              the parent node
	 * @param idStringGenerator generator for the artificial tax id, or {@code null} to skip creation
	 * @return the data child node, or {@code null} if none exists and no generator was given
	 */
	public TaxIdNode dataNode(TaxIdNode node, IDStringGenerator idStringGenerator) {
		TaxIdNode substitute = node.getDataChild();
		if (substitute == null) {
			synchronized (node) {
				substitute = node.getDataChild();
				if (substitute == null && idStringGenerator != null) {
					substitute = taxIdNodeTrie.get(idStringGenerator.generateID(nextArtCounter.getAndIncrement()), true);
					substitute.rank = (short) Rank.DATA.ordinal();
					substitute.name = "Data for " + node.getTaxId();
					node.addSubNode(substitute);
				}
			}
		}
		return substitute;
	}

	/**
	 * Returns the artificial {@code FILE}-rank child of the given node with the given
	 * name, creating it (with a freshly generated artificial tax id) if it does not yet
	 * exist. Thread-safe.
	 *
	 * @param node              the parent node
	 * @param name              the name of the file child
	 * @param idStringGenerator generator for the artificial tax id, or {@code null} to skip creation
	 * @return the file child node, or {@code null} if none exists and no generator was given
	 */
	public TaxIdNode fileNode(TaxIdNode node, String name, IDStringGenerator idStringGenerator) {
		TaxIdNode substitute = node.getChildWithName(name);
		if (substitute == null) {
			synchronized (node) {
				substitute = node.getChildWithName(name);
				if (substitute == null && idStringGenerator != null) {
					substitute = taxIdNodeTrie.get(idStringGenerator.generateID(nextArtCounter.getAndIncrement()), true);
					substitute.rank = (short) Rank.FILE.ordinal();
					substitute.name = name;
					node.addSubNode(substitute);
				}
			}
		}
		return substitute;
	}

	/**
	 * Returns the artificial {@code ID}-rank child of the given node whose name is the
	 * given byte sub-array range, creating it (with a freshly generated artificial tax
	 * id) if it does not yet exist. Thread-safe.
	 *
	 * @param node              the parent node
	 * @param array             the byte array holding the child name
	 * @param start             the start index of the name within the array (inclusive)
	 * @param end               the end index of the name within the array (exclusive)
	 * @param idStringGenerator generator for the artificial tax id
	 * @return the id child node
	 */
	public TaxIdNode idNode(TaxIdNode node, byte[] array, int start, int end, IDStringGenerator idStringGenerator) {
		TaxIdNode substitute = node.getChildWithName(array, start, end);
		if (substitute == null) {
			synchronized (node) {
				substitute = node.getChildWithName(array, start, end);
				if (substitute == null) {
					substitute = taxIdNodeTrie.get(idStringGenerator.generateID(nextArtCounter.getAndIncrement()), true);
					substitute.rank = (short) Rank.ID.ordinal();
					substitute.name = new String(array, start, end - start, StandardCharsets.UTF_8);
					node.addSubNode(substitute);
				}
			}
		}
		return substitute;
	}

	/**
	 * Recomputes the {@code position} of every node by a fresh depth-first traversal.
	 */
	public void reinitPositions() {
		root.initPositions(0, 0);
	}

	/**
	 * Sorts the given list of nodes into taxonomy-tree order in place and returns it.
	 *
	 * @param taxids the list of nodes to sort
	 * @return the same list, sorted in taxonomy-tree order
	 */
	public static List<TaxIdNode> sortNodes(List<TaxIdNode> taxids) {
		Collections.sort(taxids, NODE_COMPARATOR);
		return taxids;
	}

	/**
	 * Returns the root node of the tree.
	 *
	 * @return the root node
	 */
	protected TaxIdNode getRoot() {
		return root;
	}

	/**
	 * Returns the node with the given tax id, or {@code null} if there is none.
	 *
	 * @param taxId the tax id to look up
	 * @return the matching node, or {@code null}
	 */
	public TaxIdNode getNodeByTaxId(String taxId) {
		return taxIdNodeTrie.get(taxId);
	}

	/**
	 * Returns the node whose tax id equals the given byte sub-array range, or
	 * {@code null} if there is none.
	 *
	 * @param seq   the byte array holding the tax id
	 * @param start the start index (inclusive)
	 * @param end   the end index (exclusive)
	 * @return the matching node, or {@code null}
	 */
	public TaxIdNode getNodeByTaxId(byte[] seq, int start, int end) {
		return taxIdNodeTrie.get(seq, start, end);
	}

	/**
	 * A node of the {@link TaxTree}, holding its parent, child nodes and bookkeeping
	 * flags used while building filter lists and databases.
	 */
	public static class TaxIdNode extends TaxIdInfo {
		private static final long serialVersionUID = 1L;

		/** The direct child nodes, or {@code null} if this node has none. */
		private List<TaxIdNode> subNodes;
		/** The parent node, or {@code null} for the root. */
		private TaxIdNode parent;
		/** The depth of this node (root = 0), assigned by {@link #initPositions(int, int)}. */
		private int depth;
		/** Whether this node is required and must be retained in a {@link SmallTaxTree}. */
		private boolean required;
		/** The number of RefSeq regions associated with this node and its descendants. */
		private int refSeqRegions;

		/**
		 * Creates a node with the given tax id and no rank.
		 *
		 * @param taxId the tax id of the node
		 */
		public TaxIdNode(String taxId) {
			this(taxId, null);
		}

		/**
		 * Creates a node with the given tax id and rank.
		 *
		 * @param taxId the tax id of the node
		 * @param rank  the rank of the node, or {@code null} if unknown
		 */
		public TaxIdNode(String taxId, Rank rank) {
			super(taxId, rank);
			subNodes = null;
		}

		/**
		 * Returns the direct child with the artificial {@code DATA} rank, or {@code null}
		 * if there is none.
		 *
		 * @return the data child, or {@code null}
		 */
		public TaxIdNode getDataChild() {
			if (subNodes == null) {
				return null;
			}
			for (int i = 0; i < subNodes.size(); i++) {
				TaxIdNode child = subNodes.get(i);
				if (Rank.DATA.ordinal() == child.getRankOrdinal()) {
					return child;
				}
			}
			return null;
		}

		/**
		 * Returns the direct child whose name equals the given string, or {@code null} if none.
		 *
		 * @param name the name to match
		 * @return the matching child, or {@code null}
		 */
		public TaxIdNode getChildWithName(String name) {
			if (subNodes == null) {
				return null;
			}
			for (int i = 0; i < subNodes.size(); i++) {
				TaxIdNode child = subNodes.get(i);
				if (child.name.equals(name)) {
					return child;
				}
			}
			return null;
		}

		/**
		 * Returns the direct child whose name equals the given byte sub-array range, or
		 * {@code null} if none.
		 *
		 * @param array the byte array holding the name
		 * @param start the start index (inclusive)
		 * @param end   the end index (exclusive)
		 * @return the matching child, or {@code null}
		 */
		public TaxIdNode getChildWithName(byte[] array, int start, int end) {
			if (subNodes == null) {
				return null;
			}
			for (int i = 0; i < subNodes.size(); i++) {
				TaxIdNode child = subNodes.get(i);
				if (ByteArrayUtil.equals(array, start, end, child.name)) {
					return child;
				}
			}
			return null;
		}

		/**
		 * Assigns {@code position} values (depth-first order, starting from {@code counter}) and
		 * {@code depth} values (root = 0) to this node and its descendants, and returns the last
		 * position used.
		 *
		 * @param counter the position to assign to this node
		 * @param depth the depth to assign to this node (its distance from the root)
		 * @return the last position assigned in this subtree
		 */
		protected int initPositions(int counter, int depth) {
			position = counter;
			this.depth = depth;
			if (subNodes != null) {
				// Not using iterator for more efficiency (less object burn)
				int size = subNodes.size();
				for (int i = 0; i < size; i++) {
					counter = subNodes.get(i).initPositions(counter + 1, depth + 1);
				}
			}
			return counter;
		}

		/**
		 * Returns the depth of this node (root = 0), as assigned by the last
		 * {@link #initPositions(int, int)} traversal.
		 *
		 * @return the depth of this node
		 */
		public int getDepth() {
			return depth;
		}

		/**
		 * Increments the RefSeq region count of this node and all of its ancestors.
		 */
		public void incRefSeqRegions() {
			for (TaxIdNode node = this; node != null; node = node.parent) {
				node.refSeqRegions++;
			}
		}

		/**
		 * Returns the number of RefSeq regions associated with this node.
		 *
		 * @return the RefSeq region count
		 */
		public int getRefSeqRegions() {
			return refSeqRegions;
		}

		/**
		 * Marks this node and all of its ancestors as required, so they are retained
		 * when building a {@link SmallTaxTree}.
		 */
		public void markRequired() {
			if (required) {
				return;
			}
			required = true;
			if (parent != null) {
				parent.markRequired();
			}
		}
		
		/**
		 * Returns whether this node has been marked as required.
		 *
		 * @return {@code true} if this node is required
		 */
		public boolean isRequired() {
			return required;
		}

		/**
		 * Returns the direct child nodes, or {@code null} if there are none.
		 *
		 * @return the list of child nodes, or {@code null}
		 */
		protected List<TaxIdNode> getSubNodes() {
			return subNodes;
		}
		
		private void addSubNode(TaxIdNode node) {
			if (subNodes == null) {
				subNodes = new ArrayList<>();
			}
			node.parent = this;
			// Keep the child's depth consistent with its parent. During the initial build the parent's
			// depth is not final yet, but initPositions() overwrites every depth afterwards; for
			// artificial nodes added after that traversal (data/file/id nodes) this keeps depth correct.
			node.depth = depth + 1;
			subNodes.add(node);
		}

		public TaxIdNode getParent() {
			return parent;
		}
	}

	/**
	 * A {@link DigitTrie} mapping tax id strings to their {@link TaxIdNode}s, creating
	 * missing nodes on demand.
	 */
	public static class TaxIdNodeTrie extends DigitTrie<TaxIdNode> {
		/**
		 * Creates an empty trie.
		 */
		public TaxIdNodeTrie() {
		}

		@Override
		protected TaxIdNode createInGet(byte[] seq, int start, int end, Object createContext) {
			return new TaxIdNode(new String(seq, start, end - start, StandardCharsets.UTF_8), null);
		}

		@Override
		protected TaxIdNode createInGet(String digits, Object createContext) {
			return new TaxIdNode(digits, null);
		}
	}

	/**
	 * Generates the tax id string for an artificial node from a running counter value.
	 */
	public interface IDStringGenerator {
		/**
		 * Generates an artificial tax id string from the given counter value.
		 *
		 * @param counter the running counter value
		 * @return the generated tax id string
		 */
		public String generateID(int counter);
	}
}
