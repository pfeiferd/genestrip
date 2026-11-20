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
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;

import org.metagene.genestrip.refseq.GenomeKeyTrie;
import org.metagene.genestrip.util.ByteArrayUtil;
import org.metagene.genestrip.util.DigitTrie;

/**
 * The taxonomy tree, read from files in a {@link TaxTreeFormat}, by default the NCBI taxonomy dump
 * ({@link NCBITaxTreeFormat}). Nodes are looked up by tax id via an internal trie. The tree supports
 * ancestor and lowest-common-ancestor queries and can create artificial child nodes
 * (data, file and id nodes) used to track where k-mers originate from. It can be
 * converted into a compact {@link SmallTaxTree} for use during matching.
 */
public class TaxTree {
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

	/** File name of the NCBI taxonomy nodes dump; see {@link NCBITaxTreeFormat#NODES_DMP}. */
	public static final String NODES_DMP = NCBITaxTreeFormat.NODES_DMP;
	/** File name of the NCBI taxonomy names dump; see {@link NCBITaxTreeFormat#NAMES_DMP}. */
	public static final String NAMES_DMP = NCBITaxTreeFormat.NAMES_DMP;

	private final TaxIdNode root;
	private final TaxIdNodeTrie taxIdNodeTrie;

	/**
	 * Loads the taxonomy tree from the NCBI taxonomy dump ({@code nodes.dmp} and {@code names.dmp})
	 * located in the given directory.
	 *
	 * @param path        directory containing the NCBI taxonomy dump files
	 * @param progressBar whether to show a progress bar while reading the files
	 */
	public TaxTree(File path, boolean progressBar) {
		this(path, progressBar, new NCBITaxTreeFormat());
	}

	/**
	 * Loads the taxonomy tree from the files in the given directory, read in the given format.
	 *
	 * @param path        directory containing the taxonomy files
	 * @param progressBar whether to show a progress bar while reading the files
	 * @param format      the format of the taxonomy files
	 */
	public TaxTree(File path, boolean progressBar, TaxTreeFormat format) {
		try {
			taxIdNodeTrie = new TaxIdNodeTrie();
			root = format.read(path, new Builder(), progressBar);
			if (root == null) {
				throw new IllegalStateException("No root node in the taxonomy files of " + path);
			}

			root.initPositions(0, 0);
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}
	
	/**
	 * Returns the name an {@code OTHER} node hanging under the given taxon carries: the taxon's own
	 * tax id followed by the name of the {@link Rank#OTHER} rank.
	 * <p>
	 * It follows the convention the refinement uses for its own {@code OTHER} slot, where a concrete
	 * child contributes its tax id to a node's name and the leftover bucket contributes the word
	 * (see {@code UpdateStoreGoal.buildFirstLastName}). Both take the word from
	 * {@link Rank#getName()}, so there is one place it is written down.
	 *
	 * @param taxId the tax id of the taxon the node hangs under
	 * @return the name for its {@code OTHER} node
	 */
	public static String otherNodeName(String taxId) {
		return taxId + Rank.OTHER.getName();
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
	 * Returns the artificial {@code DATA}-rank child of the given node, creating it
	 * (with a freshly generated artificial tax id) if it does not yet exist. Thread-safe.
	 *
	 * @param node              the parent node
	 * @param idStringGenerator generator for the artificial tax id, or {@code null} to skip creation
	 * @return the data child node, or {@code null} if none exists and no generator was given
	 */
	/**
	 * Returns the artificial {@code OTHER}-rank child of the given node, creating it (with a
	 * freshly generated artificial tax id) if it does not yet exist. Thread-safe.
	 * <p>
	 * The node stands for the genomes below {@code node} that nothing in the database names. It is
	 * created empty and never receives a k-mer during the entry phase, which is what distinguishes
	 * it from the other artificial nodes.
	 *
	 * @param node              the parent node
	 * @param idStringGenerator generator for the artificial tax id, or {@code null} to skip creation
	 * @return the other child node, or {@code null} if none exists and no generator was given
	 */
	public TaxIdNode otherNode(TaxIdNode node, IDStringGenerator idStringGenerator) {
		TaxIdNode substitute = node.getOtherChild();
		if (substitute == null) {
			synchronized (node) {
				substitute = node.getOtherChild();
				if (substitute == null && idStringGenerator != null) {
					substitute = taxIdNodeTrie.get(idStringGenerator.generateID(nextArtCounter.getAndIncrement()), true);
					substitute.rank = (short) Rank.OTHER.ordinal();
					// The tax id itself still comes from the counter: these are kept in a DigitTrie,
					// which takes digits and nothing else, so it cannot carry the parent's id and a
					// word.
					substitute.name = otherNodeName(node.getTaxId());
					node.addSubNode(substitute);
				}
			}
		}
		return substitute;
	}

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
	 * Returns the artificial {@code GENOME}-rank child of the given node standing for the genome
	 * of the accession {@code array[start, end)}, creating it (with a freshly generated artificial
	 * tax id) if it does not yet exist. Thread-safe.
	 * <p>
	 * The genome is recognised from the accession alone, by
	 * {@link GenomeKeyTrie#genomeKeyLength(byte[], int, int)} -- the same rule
	 * {@code maxGenomesPerTaxid} counts genomes with, so a node is created per genome exactly where
	 * a genome is counted. Every contig of one draft assembly therefore resolves to the same child,
	 * whether or not the contigs arrived in the same fasta file, which is what distinguishes this
	 * from {@link #fileNode(TaxIdNode, String, IDStringGenerator)}.
	 *
	 * @param node              the parent node
	 * @param array             the byte array holding the accession
	 * @param start             the start index of the accession within the array (inclusive)
	 * @param end               the end index of the accession within the array (exclusive)
	 * @param idStringGenerator generator for the artificial tax id
	 * @return the genome child node
	 */
	public TaxIdNode genomeNode(TaxIdNode node, byte[] array, int start, int end, IDStringGenerator idStringGenerator) {
		int keyEnd = start + GenomeKeyTrie.genomeKeyLength(array, start, end);
		TaxIdNode substitute = node.getChildWithName(array, start, keyEnd);
		if (substitute == null) {
			synchronized (node) {
				substitute = node.getChildWithName(array, start, keyEnd);
				if (substitute == null) {
					substitute = taxIdNodeTrie.get(idStringGenerator.generateID(nextArtCounter.getAndIncrement()), true);
					substitute.rank = (short) Rank.GENOME.ordinal();
					substitute.name = new String(array, start, keyEnd - start, StandardCharsets.UTF_8);
					node.addSubNode(substitute);
				}
			}
		}
		return substitute;
	}

	/**
	 * Returns the artificial {@code GENOME}-rank child of the given node with the given name, creating
	 * it (with a freshly generated artificial tax id) if it does not yet exist. Thread-safe.
	 * <p>
	 * For a genome the accession does not identify, which is then named by its taxon instead; see
	 * {@code RefSeqCategory#isOneGenomePerTaxon()}.
	 *
	 * @param node              the parent node
	 * @param name              the name of the genome
	 * @param idStringGenerator generator for the artificial tax id
	 * @return the genome child node
	 */
	public TaxIdNode genomeNode(TaxIdNode node, String name, IDStringGenerator idStringGenerator) {
		TaxIdNode substitute = node.getChildWithName(name);
		if (substitute == null) {
			synchronized (node) {
				substitute = node.getChildWithName(name);
				if (substitute == null) {
					substitute = taxIdNodeTrie.get(idStringGenerator.generateID(nextArtCounter.getAndIncrement()), true);
					substitute.rank = (short) Rank.GENOME.ordinal();
					substitute.name = name;
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
		/** The number of RefSeq contigs associated with this node and its descendants. */
		private int refSeqContigs;
		/** Whether the catalog files this taxon's sequences under a category with one genome per taxon. */
		private boolean oneGenomePerTaxon;

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
		/**
		 * Returns the direct {@code OTHER}-rank child of this node, or {@code null} if none.
		 *
		 * @return the other child, or {@code null}
		 */

		public TaxIdNode getOtherChild() {
			if (subNodes == null) {
				return null;
			}
			for (int i = 0; i < subNodes.size(); i++) {
				TaxIdNode child = subNodes.get(i);
				if (Rank.OTHER.ordinal() == child.getRankOrdinal()) {
					return child;
				}
			}
			return null;
		}

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
		 * Increments the RefSeq contig count of this node and all of its ancestors.
		 */
		public void incRefSeqContigs() {
			for (TaxIdNode node = this; node != null; node = node.parent) {
				node.refSeqContigs++;
			}
		}

		/**
		 * Marks this taxon as one whose sequences the catalog files under a category with one genome
		 * per taxon; see {@code RefSeqCategory#isOneGenomePerTaxonRow(byte[], int, int)}.
		 */
		public void markOneGenomePerTaxon() {
			oneGenomePerTaxon = true;
		}

		/**
		 * Returns whether this taxon is marked by {@link #markOneGenomePerTaxon()}.
		 *
		 * @return whether this taxon has one genome per taxon
		 */
		public boolean isOneGenomePerTaxon() {
			return oneGenomePerTaxon;
		}

		/**
		 * Returns the number of RefSeq contigs associated with this node.
		 *
		 * @return the RefSeq contig count
		 */
		public int getRefSeqContigs() {
			return refSeqContigs;
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
	 * Builds the nodes of a tree while a {@link TaxTreeFormat} reads its files, so that a format
	 * needs no access to the internals of {@link TaxIdNode}. Only the tree hands out instances.
	 */
	public final class Builder {
		private Builder() {
		}

		/**
		 * Returns the node with the tax id given by the byte sub-array range, optionally creating it.
		 *
		 * @param seq    the byte array holding the tax id
		 * @param start  the start index (inclusive)
		 * @param end    the end index (exclusive)
		 * @param create whether to create the node if it does not exist yet
		 * @return the node, or {@code null} if it does not exist and {@code create} is {@code false}
		 */
		public TaxIdNode getNode(byte[] seq, int start, int end, boolean create) {
			// Any non-null context makes the trie create the node, Boolean.FALSE included.
			return taxIdNodeTrie.get(seq, start, end, create ? Boolean.TRUE : null);
		}

		/**
		 * Returns the node with the given tax id, optionally creating it.
		 *
		 * @param taxId  the tax id
		 * @param create whether to create the node if it does not exist yet
		 * @return the node, or {@code null} if it does not exist and {@code create} is {@code false}
		 */
		public TaxIdNode getNode(String taxId, boolean create) {
			return taxIdNodeTrie.get(taxId, create ? Boolean.TRUE : null);
		}

		/**
		 * Makes {@code child} a child of {@code parent}.
		 *
		 * @param parent the parent node
		 * @param child  the child node
		 */
		public void addSubNode(TaxIdNode parent, TaxIdNode child) {
			parent.addSubNode(child);
		}

		/**
		 * Sets the rank of the node.
		 *
		 * @param node the node
		 * @param rank the rank, or {@code null} if unknown
		 */
		public void setRank(TaxIdNode node, Rank rank) {
			node.rank = (short) (rank == null ? -1 : rank.ordinal());
		}

		/**
		 * Sets the name of the node.
		 *
		 * @param node the node
		 * @param name the name
		 */
		public void setName(TaxIdNode node, String name) {
			node.name = name;
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
