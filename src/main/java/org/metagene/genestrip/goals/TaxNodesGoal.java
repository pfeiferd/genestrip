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
package org.metagene.genestrip.goals;

import java.io.IOException;
import java.util.HashSet;
import java.util.Set;

import org.metagene.genestrip.GSConfigKey;
import org.metagene.genestrip.GSGoalKey;
import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.tax.Rank;
import org.metagene.genestrip.tax.TaxIdCollector;
import org.metagene.genestrip.tax.TaxTree;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;

public class TaxNodesGoal extends ObjectGoal<Set<TaxIdNode>, GSProject> {
	private final ObjectGoal<TaxTree, GSProject> taxTreeGoal;

	@SafeVarargs
	public TaxNodesGoal(GSProject project, ObjectGoal<TaxTree, GSProject> taxTreeGoal,
			Goal<GSProject>... deps) {
		super(project, GSGoalKey.TAXNODES, Goal.append(deps, taxTreeGoal));
		this.taxTreeGoal = taxTreeGoal;
	}
	
	public ObjectGoal<TaxTree, GSProject> getTaxTreeGoal() {
		return taxTreeGoal;
	}

	@Override
	protected void doMakeThis() {
		try {
			TaxIdCollector taxIdCollector = new TaxIdCollector(taxTreeGoal.get());
			Set<TaxIdNode> excludes = new HashSet<TaxIdNode>();
			Set<TaxIdNode> taxIdNodes = taxIdCollector.readFromFile(getProject().getTaxIdsFile(), excludes);
			if (getLogger().isDebugEnabled()) {
				getLogger().debug("Requested tax ids: " + taxIdNodes);
				getLogger().debug("Excluded tax ids: " + excludes);
			}
			taxIdNodes = taxIdCollector.withDescendants(taxIdNodes, (Rank) configValue(GSConfigKey.RANK_COMPLETION_DEPTH));
			excludes = taxIdCollector.withDescendants(excludes, null);
			taxIdNodes.removeAll(excludes);
			if (getLogger().isDebugEnabled()) {
				getLogger().debug("Number of completed tax ids: " + taxIdNodes.size());
			}
			
			set(taxIdNodes);
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}
}
