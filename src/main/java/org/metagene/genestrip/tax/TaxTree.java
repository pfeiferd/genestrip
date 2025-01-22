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
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import org.metagene.genestrip.io.BufferedLineReader;
import org.metagene.genestrip.io.StreamProvider;
import org.metagene.genestrip.util.ByteArrayUtil;
import org.metagene.genestrip.util.DigitTrie;

public class TaxTree {
	private static int MAX_LINE_SIZE = 4096;

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

	public static final String NODES_DMP = "nodes.dmp";
	public static final String NAMES_DMP = "names.dmp";

	private final TaxIdNode root;
	private final TaxIdNodeTrie taxIdNodeTrie;

	public TaxTree(File path) {
		try {
			taxIdNodeTrie = new TaxIdNodeTrie();
			try (InputStream is = createNodesResource(path)) {
				root = readNodesFromStream(is);
			}
			try (InputStream is = createNamesResource(path)) {
				readNamesFromStream(is);
			}
			initPositions(0, root);
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}
	
	public SmallTaxTree toSmallTaxTree() {
		return new SmallTaxTree(this);
	}

	public boolean isAncestorOf(TaxIdNode node, TaxIdNode ancestor) {
		while (node != null) {
			if (node.equals(ancestor)) {
				return true;
			}
			node = node.parent;
		}

		return false;
	}

	public TaxIdNode getLeastCommonAncestor(final TaxIdNode node1, final TaxIdNode node2) {
		for (TaxIdNode res = node1; res != null; res = res.parent) {
			for (TaxIdNode ancestor2 = node2; ancestor2 != null; ancestor2 = ancestor2.parent) {
				if (res == ancestor2) {
					return res;
				}
			}
		}
		return null;
	}

	protected InputStream createNodesResource(File path) throws IOException {
		return StreamProvider.getInputStreamForFile(new File(path, NODES_DMP));
	}

	protected InputStream createNamesResource(File path) throws IOException {
		return StreamProvider.getInputStreamForFile(new File(path, NAMES_DMP));
	}

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
						String name = new String(target, a + 2, b - a - 3);
						if (node.name == null || scientific) {
							node.name = name;
						}
					}
				}
			}
		}
	}

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

	protected int initPositions(int counter, TaxIdNode taxIdNode) {
		taxIdNode.position = counter;
		List<TaxIdNode> subNodes = taxIdNode.subNodes;
		if (subNodes != null) {
			// Not using iterator for more efficiency (less object burn)
			int size = subNodes.size();
			for (int i = 0; i < size; i++) {
				counter = initPositions(counter + 1, subNodes.get(i));
			}
		}
		return counter;
	}

	public static List<TaxIdNode> sortNodes(List<TaxIdNode> taxids) {
		Collections.sort(taxids, NODE_COMPARATOR);
		return taxids;
	}

	protected TaxIdNode getRoot() {
		return root;
	}

	public TaxIdNode getNodeByTaxId(String taxId) {
		return taxIdNodeTrie.get(taxId);
	}

	public TaxIdNode getNodeByTaxId(byte[] seq, int start, int end) {
		return taxIdNodeTrie.get(seq, start, end);
	}

	public static class TaxIdNode extends TaxIdInfo {
		private static final long serialVersionUID = 1L;

		private List<TaxIdNode> subNodes;
		private TaxIdNode parent;
		private boolean required;
		private int refSeqRegions;

		public TaxIdNode(String taxId) {
			this(taxId, null);
		}

		public TaxIdNode(String taxId, Rank rank) {
			super(taxId, rank);
			subNodes = null;
		}

		public void incRefSeqRegions() {
			for (TaxIdNode node = this; node != null; node = node.parent) {
				node.refSeqRegions++;
			}
		}

		public int getRefSeqRegions() {
			return refSeqRegions;
		}

		public void markRequired() {
			if (required) {
				return;
			}
			required = true;
			if (parent != null) {
				parent.markRequired();
			}
		}
		
		public boolean isRequired() {
			return required;
		}

		protected List<TaxIdNode> getSubNodes() {
			return subNodes;
		}
		
		private void addSubNode(TaxIdNode node) {
			if (subNodes == null) {
				subNodes = new ArrayList<TaxIdNode>();
			}
			node.parent = this;
			subNodes.add(node);
		}

		public TaxIdNode getParent() {
			return parent;
		}
	}

	public static class TaxIdNodeTrie extends DigitTrie<TaxIdNode> {
		@Override
		protected TaxIdNode createInGet(byte[] seq, int start, int end, Object createContext) {
			return new TaxIdNode(new String(seq, start, end - start), null);
		}
	}
}
