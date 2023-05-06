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
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.metagene.genestrip.util.StreamProvider;

public class TaxTree {
	public enum Rank {
		NO_RANK("no rank"), SUPERKINGDOM("superkingdom"), KINGDOM("kingdom"), PHYLUM("phylum"), SUBPHYLUM("subphylum"),
		SUPERCLASS("superclass"), CLASS("class"), SUBCLASS("subclass"), SUPERORDER("superorder"), ORDER("order"),
		SUBORDER("suborder"), SUPERFAMILY("superfamily"), FAMILY("family"), SUBFAMILY("subfamily"), CLADE("clade"),
		GENUS("genus"), SUBGENUS("subgenus"), SPECIES_GROUP("species group"), SPECIES("species"), VARIETAS("varietas"),
		SUBSPECIES("subspecies"), STRAIN("strain");

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
	private final Map<String, TaxIdNode> taxIdToNode;

	public TaxTree(File path) {
		try {
			taxIdToNode = new HashMap<String, TaxIdNode>();
			root = readNodesFromStream(createNodesResource(path), taxIdToNode);
			readNamesFromStream(createNamesResource(path), taxIdToNode);
			initPositions(0, root);
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
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

	protected InputStream createNodesResource(File path) throws IOException {
		return StreamProvider.getInputStreamForFile(new File(path, NODES_DMP));
	}

	protected InputStream createNamesResource(File path) throws IOException {
		return StreamProvider.getInputStreamForFile(new File(path, NAMES_DMP));
	}

	protected void readNamesFromStream(InputStream stream, Map<String, TaxIdNode> taxIdToNode) throws IOException {
		try (BufferedReader br = new BufferedReader(new InputStreamReader(stream, StandardCharsets.UTF_8))) {
			String line = null;
			while ((line = br.readLine()) != null) {
				boolean scientific = line.indexOf("scientific name") != -1;
				int a = line.indexOf("|");
				int b = line.indexOf("|", a + 1);
				if (a != -1 && b != -1) {
					String taxA = line.substring(0, a).trim();
					String name = line.substring(a + 1, b).trim();

					TaxIdNode node = taxIdToNode.get(taxA);
					if (node != null) {
						if (node.name == null || scientific) {
							node.name = name;
						}
					}
				}
			}
		}
	}

	protected TaxIdNode readNodesFromStream(InputStream stream, Map<String, TaxIdNode> taxIdToNode) throws IOException {
		TaxIdNode res = null;
		try (BufferedReader br = new BufferedReader(new InputStreamReader(stream, StandardCharsets.UTF_8))) {
			String line = null;
			while ((line = br.readLine()) != null) {
				int a = line.indexOf("|");
				int b = line.indexOf("|", a + 1);
				int c = line.indexOf("|", b + 1);
				if (a != -1 && b != -1) {
					String taxA = line.substring(0, a).trim();
					String taxB = line.substring(a + 1, b).trim();
					String level = line.substring(b + 1, c).trim();

					TaxIdNode nodeA;
					if (taxA.equals(taxB) && taxA.equals("1")) {
						res = nodeA = new TaxIdNode(taxA, null, Rank.byName(level));
					} else {
						nodeA = new TaxIdNode(taxA, taxB, Rank.byName(level));
					}
					taxIdToNode.put(taxA, nodeA);
				}
			}
		}

		for (TaxIdNode node : taxIdToNode.values()) {
			TaxIdNode parent = taxIdToNode.get(node.getParentTaxId());
			if (parent != null) {
				parent.addSubNode(node);
				node.parent = parent;
			}
		}

		return res;
	}

	protected int initPositions(int counter, TaxIdNode taxIdNode) {
		taxIdNode.position = counter;
		for (TaxIdNode subNode : taxIdNode.subNodes) {
			counter = initPositions(counter + 1, subNode);
		}
		return counter;
	}

	public TaxIdNode getRoot() {
		return root;
	}

	public TaxIdNode getNodeByTaxId(String taxId) {
		return taxIdToNode.get(taxId);
	}

	public static class TaxIdNode implements Comparable<TaxIdNode> {
		private final String taxId;
		private final String parentTaxId;
		private String name;
		private final List<TaxIdNode> subNodes;
		private TaxIdNode parent;
		private int position;
		private Rank rank;

		public TaxIdNode(String taxId, String name, Rank rank) {
			this.parentTaxId = name;
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

		public String getParentTaxId() {
			return parentTaxId;
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
}
