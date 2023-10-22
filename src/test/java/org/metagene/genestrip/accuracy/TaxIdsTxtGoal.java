package org.metagene.genestrip.accuracy;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;

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

	@SafeVarargs
	public TaxIdsTxtGoal(GSProject project, String name, ObjectGoal<TaxTree, GSProject> taxTreeGoal,
			Goal<GSProject>... dependencies) {
		super(project, name, new File(project.getProjectDir(), "taxids.txt"), append(dependencies, taxTreeGoal));
		this.taxTreeGoal = taxTreeGoal;
	}

	@Override
	protected void makeFile(File file) throws IOException {
		PrintStream out = new PrintStream(StreamProvider.getOutputStreamForFile(file));

		TaxTree taxTree = taxTreeGoal.get();
		for (String taxid : FastaTransformGoal.TAXIDS) {
			TaxIdNode node = taxTree.getNodeByTaxId(taxid);
			if (node != null) {
				TaxIdNode genusNode = AccuracyMatchGoal.getRankTaxNode(taxTree, taxid, Rank.SPECIES);
				if (genusNode != null) {
					out.println(genusNode.getTaxId());
				}
				else {
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
		out.close();
	}
}
