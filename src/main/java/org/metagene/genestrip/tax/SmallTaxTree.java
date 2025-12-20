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
import org.metagene.genestrip.util.DigitTrie;

public class SmallTaxTree implements Serializable, Iterable<SmallTaxTree.SmallTaxIdNode> {
	private static final long serialVersionUID = 1L;

	private transient Comparator<String> taxIdComparator;
	private transient int countSize;
	private transient SmallTaxIdNode root;
	private transient DigitTrie<SmallTaxIdNode> taxIdNodeTrie;
	private transient Object owner;

	public SmallTaxTree(TaxTree taxTree) {
		root = taxTree.getRoot() == null ? null : new SmallTaxIdNode(taxTree.getRoot());
		taxIdNodeTrie = new DigitTrie<SmallTaxIdNode>();
		root.initTrie(taxIdNodeTrie);
	}

	private void writeObject(ObjectOutputStream out) throws IOException {
		root.writeTree(out);
	}

	private void readObject(ObjectInputStream in) throws IOException, ClassNotFoundException {
		root = SmallTaxIdNode.readTree(in);
		taxIdNodeTrie = new DigitTrie<SmallTaxIdNode>();
		root.initTrie(taxIdNodeTrie);
	}

	public SmallTaxIdNode setSubNodes(String taxId, SmallTaxIdNode[] subNodes) {
		SmallTaxIdNode node = taxIdNodeTrie.get(taxId);
		if (node != null) {
			node.setSubNodes(subNodes);
			node.initTrie(taxIdNodeTrie);
		}

		return node;
	}

	public void reinitPositions() {
		root.initPositions(0);
	}

	public void initCountSize(int countSize) {
		if (countSize <= 0) {
			throw new IllegalArgumentException("Initialization size must by >= 0.");
		}
		if (this.countSize > 0 && countSize != this.countSize) {
			throw new IllegalStateException("Count count size can only be initialized once.");
		}
		this.countSize = countSize;
	}

	public void resetCounts(Object owner) {
		if (owner == null) {
			throw new NullPointerException("owner must not be null");
		}
		if (this.owner == null) {
			this.owner = owner;
		} else if (this.owner != owner) {
			throw new IllegalArgumentException("Tax tree not relaease by previous owner: " + this.owner);
		}
		root.resetCounts();
	}

	public void releaseOwner() {
		this.owner = null;
	}

	// Made final for potential inlining by JVM
	public final void incCount(final SmallTaxIdNode node, final int index, final long initKey) {
		node.incCount(index, initKey, countSize);
	}

	// Made final for potential inlining by JVM
	public final short sumCounts(SmallTaxIdNode node, final int index, final long initKey) {
		short res = 0;
		while (node != null) {
			if (node.counts != null && node.countsInitKeys != null && node.countsInitKeys[index] == initKey) {
				res += node.counts[index];
			}
			node = node.parent;
		}
		return res;
	}

	public int getCountSize() {
		return countSize;
	}

	// Made final for potential inlining by JVM
	public final boolean isAncestorOf(SmallTaxIdNode node, final SmallTaxIdNode ancestor) {
		while (node != null) {
			if (node.equals(ancestor)) {
				return true;
			}
			node = node.parent;
		}

		return false;
	}

	// Made final for potential inlining by JVM
	public final SmallTaxIdNode getLeastCommonAncestor(final SmallTaxIdNode node1, final SmallTaxIdNode node2) {
		for (SmallTaxIdNode res = node1; res != null; res = res.parent) {
			for (SmallTaxIdNode ancestor2 = node2; ancestor2 != null; ancestor2 = ancestor2.parent) {
				if (res == ancestor2) {
					return res;
				}
			}
		}
		return null;
	}

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

	public SmallTaxIdNode getRoot() {
		return root;
	}

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


	public SmallTaxIdNode getNodeByTaxId(String taxId) {
		return taxIdNodeTrie.get(taxId);
	}

	public static class SmallTaxIdNode extends TaxIdInfo {
		private static final long serialVersionUID = 1L;

		private SmallTaxIdNode[] subNodes;
		private boolean requested;

		protected transient SmallTaxIdNode parent;
		private transient short[] counts;
		private transient long[] countsInitKeys;
		// Made public for inlining
		public transient short storeIndex;

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

		public boolean isRequested() {
			return requested;
		}
		
		public void setRequested(boolean value) {
			requested = value;
		}

		public SmallTaxIdNode getParent() {
			return parent;
		}

		public final short getStoreIndex() {
			return storeIndex;
		}

		public void setStoreIndex(short storeIndex) {
			this.storeIndex = storeIndex;
		}

		public SmallTaxIdNode[] getSubNodes() {
			return subNodes;
		}

		public void setSubNodes(SmallTaxIdNode[] subNodes) {
			this.subNodes = subNodes;
			for (SmallTaxIdNode subNode : subNodes) {
				subNode.parent = this;
			}
		}

		public final void incCount(final int index, final long initKey, final int size) {
			if (counts == null || countsInitKeys == null) {
				synchronized (this) {
					if (counts == null) {
						counts = new short[size];
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
