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

		for (String taxid : accessionNumber2TaxidGoal.get().values()) {
			addNodesForTaxid(taxid, nodes);
		}

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
