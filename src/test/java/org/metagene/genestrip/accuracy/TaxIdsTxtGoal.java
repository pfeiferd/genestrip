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
package org.metagene.genestrip.accuracy;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.io.StreamProvider;
import org.metagene.genestrip.make.FileListGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.tax.TaxTree;
import org.metagene.genestrip.tax.TaxTree.Rank;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;

public class TaxIdsTxtGoal extends FileListGoal<GSProject> {
	private final ObjectGoal<TaxTree, GSProject> taxTreeGoal;

	private final ObjectGoal<Map<String, String>, GSProject> accessionNumber2TaxidGoal;

	@SafeVarargs
	public TaxIdsTxtGoal(GSProject project, String name, ObjectGoal<TaxTree, GSProject> taxTreeGoal,
			ObjectGoal<Map<String, String>, GSProject> accessionNumber2TaxidGoal, Goal<GSProject>... dependencies) {
		super(project, name, new File(project.getProjectDir(), "taxids.txt"),
				append(dependencies, taxTreeGoal, accessionNumber2TaxidGoal));
		this.taxTreeGoal = taxTreeGoal;
		this.accessionNumber2TaxidGoal = accessionNumber2TaxidGoal;
	}

	@Override
	protected void makeFile(File file) throws IOException {

		List<TaxIdNode> nodes = new ArrayList<TaxIdNode>();

		for (String taxid : FastaTransformGoal.TAXIDS) {
			addNodesForTaxid(taxid, nodes);
		}

//		for (String taxid : accessionNumber2TaxidGoal.get().values()) {
//			addNodesForTaxid(taxid, nodes);
//		}

		PrintStream out = new PrintStream(StreamProvider.getOutputStreamForFile(file));
		TaxTree.sortNodes(nodes);
		for (TaxIdNode node : nodes) {
			out.println(node.getTaxId());
		}

		out.close();
	}

	protected void addNodesForTaxid(String taxid, List<TaxIdNode> nodes) {
		TaxTree taxTree = taxTreeGoal.get();
		TaxIdNode node = taxTree.getNodeByTaxId(taxid);
		if (node != null) {
			if (!nodes.contains(node)) {
				nodes.add(node);
			}
			TaxIdNode genusNode = taxTree.getRankedNode(taxid, Rank.SPECIES);
			if (genusNode != null) {
				if (!nodes.contains(genusNode)) {
					nodes.add(genusNode);
				}
			} else {
				if (getLogger().isWarnEnabled()) {
					getLogger().warn("Missing genus for taxid " + taxid);
				}
			}
		} else {
			if (getLogger().isWarnEnabled()) {
				getLogger().warn("Missing node for taxid " + taxid);
			}
		}
	}
}
