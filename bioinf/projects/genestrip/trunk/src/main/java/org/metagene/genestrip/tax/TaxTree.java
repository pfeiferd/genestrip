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

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.io.Serializable;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.metagene.genestrip.util.BufferedLineReader;
import org.metagene.genestrip.util.ByteArrayUtil;
import org.metagene.genestrip.util.DigitTrie;
import org.metagene.genestrip.util.StreamProvider;

public class TaxTree implements Serializable {
	private static final long serialVersionUID = 1L;

	public enum Rank {
		SUPERKINGDOM("superkingdom"), KINGDOM("kingdom"), PHYLUM("phylum"), SUBPHYLUM("subphylum"),
		SUPERCLASS("superclass"), CLASS("class"), SUBCLASS("subclass"), SUPERORDER("superorder"), ORDER("order"),
		SUBORDER("suborder"), SUPERFAMILY("superfamily"), FAMILY("family"), SUBFAMILY("subfamily"), CLADE("clade"),
		GENUS("genus"), SUBGENUS("subgenus"), SPECIES_GROUP("species group"), SPECIES("species"), VARIETAS("varietas"),
		SUBSPECIES("subspecies"), STRAIN("strain"), SEROTYPE("serotype"), FORMA("forma"),
		FORMA_SPECIALIS("forma specialis"), ISOLATE("isolate"), NO_RANK("no rank");

		private String name;

		private Rank(String name) {
			this.name = name;
		}

		public String getName() {
			return name;
		}

		public static Rank byName(String name) {
			for (Rank r : Rank.values()) {
				if (r.name.equals(name)) {
					return r;
				}
			}
			return null;
		}

		public boolean isBelowOrEqual(Rank rank) {
			return this.ordinal() >= rank.ordinal();
		}

		public boolean isBelow(Rank rank) {
			return this.ordinal() > rank.ordinal();
		}
	}

	public static final String NODES_DMP = "nodes.dmp";
	public static final String NAMES_DMP = "names.dmp";

	private final TaxIdNode root;
//	private final Map<String, TaxIdNode> taxIdToNode;
	private final TaxIdNodeTrie taxIdNodeTrie;

	public TaxTree(File path) {
		try {
			taxIdNodeTrie = new TaxIdNodeTrie();
			root = readNodesFromStream(createNodesResource(path));
			readNamesFromStream(createNamesResource(path));
			initPositions(0, root);
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}

	public TaxIdNodeTrie getTaxIdNodeTrie() {
		return taxIdNodeTrie;
	}

	public boolean isAncestorOf(TaxIdNode node, TaxIdNode ancestor) {
		while (node != null) {
			if (node.equals(ancestor)) {
				return true;
			}
			node = node.getParent();
		}

		return false;
	}

	public TaxIdNode getLeastCommonAncestor(TaxIdNode node1, TaxIdNode node2) {
		for (; node1 != null; node1 = node1.getParent()) {
			for (TaxIdNode ancestor2 = node2; ancestor2 != null; ancestor2 = ancestor2.getParent()) {
				if (node1 == ancestor2) {
					return node1;
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
		BufferedLineReader br = new BufferedLineReader(stream);
		int size;
		byte[] target = new byte[2000];

		while ((size = br.nextLine(target)) >= 0) {
			boolean scientific = ByteArrayUtil.indexOf(target, 0, size, "scientific name") != -1;
			int a = ByteArrayUtil.indexOf(target, 0, size, '|');
			int b = ByteArrayUtil.indexOf(target, a + 1, size, '|');
			if (a != -1 && b != -1) {
				TaxIdNode node = taxIdNodeTrie.get(target, 0, a);
				if (node != null) {
					String name = new String(target, a + 1, b - a - 1);
					if (node.name == null || scientific) {
						node.name = name;
					}
				}
			}
		}
		br.close();
	}

	protected TaxIdNode readNodesFromStream(InputStream stream) throws IOException {
		TaxIdNode res = null;
		int size;
		byte[] target = new byte[2000];

		BufferedLineReader br = new BufferedLineReader(stream);
		while ((size = br.nextLine(target)) >= 0) {
			int a = ByteArrayUtil.indexOf(target, 0, size, '|');
			int b = ByteArrayUtil.indexOf(target, a + 1, size, '|');
			int c = ByteArrayUtil.indexOf(target, b + 1, size, '|');
			if (a != -1 && b != -1) {
				TaxIdNode nodeA = taxIdNodeTrie.get(target, 0, a, true);
				TaxIdNode nodeB = taxIdNodeTrie.get(target, a + 1, b, true);
				
				nodeA.parent = nodeB;
				nodeB.addSubNode(nodeA);
				String level = new String(target, b + 1, c - b - 1);
				nodeA.rank = Rank.byName(level);

				if (nodeA == nodeB && "1".equals(nodeA.name)) {
					res = nodeA;
				} 
			}
		}
		br.close();

		return res;
	}

	protected int initPositions(int counter, TaxIdNode taxIdNode) {
		taxIdNode.position = counter;
		for (TaxIdNode subNode : taxIdNode.subNodes) {
			counter = initPositions(counter + 1, subNode);
		}
		return counter;
	}

	public List<String> sortTaxidsViaTree(List<String> taxids) {
		Collections.sort(taxids, new Comparator<String>() {
			@Override
			public int compare(String o1, String o2) {
				TaxIdNode a = getNodeByTaxId(o1);
				TaxIdNode b = getNodeByTaxId(o2);

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
		});
		return taxids;
	}

	public static List<TaxIdNode> sortNodes(List<TaxIdNode> taxids) {
		Collections.sort(taxids, new Comparator<TaxIdNode>() {
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
		});
		return taxids;
	}

	public TaxIdNode getRoot() {
		return root;
	}

	public TaxIdNode getNodeByTaxId(String taxId) {
		return taxIdNodeTrie.get(taxId);
	}

	public static class TaxIdNode implements Comparable<TaxIdNode>, Serializable {
		private static final long serialVersionUID = 1L;

		private final String taxId;
		private String name;
		private final List<TaxIdNode> subNodes;
		private TaxIdNode parent;
		private int position;
		private Rank rank;

		public TaxIdNode(String taxId) {
			this.taxId = taxId;
			subNodes = new ArrayList<TaxIdNode>();
		}

		public TaxIdNode(String taxId, Rank rank) {
			this.taxId = taxId;
			this.rank = rank;
			subNodes = new ArrayList<TaxIdNode>();
		}

		public Rank getRank() {
			return rank;
		}

		public TaxIdNode getParent() {
			return parent;
		}

		public int getPosition() {
			return position;
		}

		public String getName() {
			return name;
		}

		public String getTaxId() {
			return taxId;
		}

		public List<TaxIdNode> getSubNodes() {
			return Collections.unmodifiableList(subNodes);
		}

		public void addSubNode(TaxIdNode node) {
			subNodes.add(node);
		}

		public void addSubNode(int position, TaxIdNode node) {
			subNodes.add(position, node);
		}

		@Override
		public int compareTo(TaxIdNode o) {
			return position - o.position;
		}

		@Override
		public String toString() {
			return "Node: " + taxId;
		}

		public TaxIdNode shallowCopy() {
			TaxIdNode copy = new TaxIdNode(taxId, rank);
			copy.name = name;
			copy.position = position;
			return copy;
		}
	}

	public void writeSortedSpecies(PrintWriter writer) throws IOException {
		writeSortedSpecies(root, writer);
	}

	protected void writeSortedSpecies(TaxIdNode node, PrintWriter writer) throws IOException {
		writer.print(node.position);
		writer.print('|');
		writer.print(node.taxId);
		writer.print('|');
		writer.println(node.name);
		for (TaxIdNode subNode : node.subNodes) {
			writeSortedSpecies(subNode, writer);
		}
	}

	public static void main(String[] args) throws IOException {
		File file = new File("sortedSpecies.txt");
		new TaxTree(new File("db/")).writeSortedSpecies(new PrintWriter(file));
	}

	public static class TaxIdNodeTrie extends DigitTrie<TaxIdNode> {
		private static final long serialVersionUID = 1L;

		@Override
		protected TaxIdNode createInGet(byte[] seq, int start, int end) {
			return new TaxIdNode(new String(seq, start, end - start));
		}
	}
}
