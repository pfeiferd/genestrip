package org.metagene.genestrip.goals;

import java.io.File;
import java.io.IOException;
import java.util.Set;

import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.tax.TaxIdCollector;
import org.metagene.genestrip.tax.TaxTree;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;
import org.metagene.genestrip.util.ArraysUtil;

public class TaxNodesGoal extends ObjectGoal<Set<TaxIdNode>, GSProject> {
	private final ObjectGoal<TaxTree, GSProject> taxTreeGoal;

	@SafeVarargs
	public TaxNodesGoal(GSProject project, String name, ObjectGoal<TaxTree, GSProject> taxTreeGoal,
			Goal<GSProject>... deps) {
		super(project, name, ArraysUtil.append(deps, taxTreeGoal));
		this.taxTreeGoal = taxTreeGoal;
	}

	@Override
	public void makeThis() {
		try {
			TaxIdCollector taxIdCollector = new TaxIdCollector(taxTreeGoal.get());
			Set<TaxIdNode> taxIdNodes = taxIdCollector.readFromFile(getProject().getTaxIdsFile());
			if (getLogger().isInfoEnabled()) {
				getLogger().info("Requested tax ids: " + taxIdNodes);
			}
			taxIdNodes = taxIdCollector.withDescendants(taxIdNodes);
			if (getLogger().isInfoEnabled()) {
//				getLogger().info("Completed tax ids: " + taxIdNodes);
				getLogger().info("Number of completed tax ids: " + taxIdNodes.size());
			}

			File filterFile = getProject().getTaxIdsFilterFile();
			if (filterFile.exists()) {
				Set<TaxIdNode> filterNodes = taxIdCollector.readFromFile(filterFile);
				taxIdNodes = taxIdCollector.restrictToAncestors(filterNodes, taxIdNodes);
				if (getLogger().isInfoEnabled()) {
//					getLogger().info("Filtered tax ids: " + taxIdNodes);
					getLogger().info("Number of completed and filtered tax ids: " + taxIdNodes.size());
				}
			}
			set(taxIdNodes);
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}

}