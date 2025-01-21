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
import java.nio.charset.StandardCharsets;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.metagene.genestrip.io.StreamProvider;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;

public class TaxIdCollector {
	private final TaxTree taxTree;

	public TaxIdCollector(TaxTree taxTree) {
		this.taxTree = taxTree;
	}

	public Set<TaxIdNode> readFromFile(File file) throws IOException {
		Set<TaxIdNode> res = new HashSet<TaxIdNode>();

		try (InputStream stream = StreamProvider.getInputStreamForFile(file);
				BufferedReader br = new BufferedReader(new InputStreamReader(stream, StandardCharsets.UTF_8))) {
			String line = null;
			while ((line = br.readLine()) != null) {
				line = line.trim();
				// Ignore whole line comments.
				if (!line.startsWith("#")) {
					// Ignore end of line comments.
					int commentIndex = line.indexOf('#');
					if (commentIndex != -1) {
						line = line.substring(0, commentIndex);
					}
					int tab = line.lastIndexOf('\t');
					String taxId;
					if (tab != -1) {
						taxId = line.substring(tab, line.length()).trim();
					} else {
						taxId = line.trim();
					}
					if (!taxId.isEmpty()) {
						TaxIdNode node = taxTree.getNodeByTaxId(taxId);
						if (node != null) {
							res.add(node);
						}
					}
				}
			}
		}

		return res;
	}

	public Set<TaxIdNode> withDescendants(Set<TaxIdNode> taxIds, Rank depth) {
		Set<TaxIdNode> res = new HashSet<TaxTree.TaxIdNode>();

		for (TaxIdNode node : taxIds) {
			completeFilterlist(res, node, depth);
		}
		return res;
	}

	public Set<TaxIdNode> restrictToAncestor(TaxIdNode ancestor, Set<TaxIdNode> taxIds) {
		Set<TaxIdNode> res = new HashSet<TaxTree.TaxIdNode>();

		for (TaxIdNode node : taxIds) {
			if (taxTree.isAncestorOf(node, ancestor)) {
				res.add(node);
			}
		}
		return res;
	}

	public Set<TaxIdNode> restrictToAncestors(Set<TaxIdNode> ancestors, Set<TaxIdNode> taxIds) {
		Set<TaxIdNode> res = new HashSet<TaxTree.TaxIdNode>();

		for (TaxIdNode node : taxIds) {
			for (TaxIdNode ancestor : ancestors) {
				if (taxTree.isAncestorOf(node, ancestor)) {
					res.add(node);
				}
			}
		}
		return res;
	}

//	public static List<TaxIdNode> sortNodes(Set<TaxIdNode> nodes) {
//		return TaxTree.sortNodes(new ArrayList<TaxIdNode>(nodes));
//	}

	private void completeFilterlist(Set<TaxIdNode> filter, TaxIdNode node, Rank depth) {
		if (node != null) {
			filter.add(node);
			List<TaxIdNode> subNodes = node.getSubNodes();
			if (subNodes != null) {
				for (TaxIdNode subNode : subNodes) {
					if (depth == null
							|| (subNode != null && subNode.getRank() != null && !subNode.getRank().isBelow(depth))) {
						completeFilterlist(filter, subNode, depth);
					}
				}
			}
		}
	}
}
