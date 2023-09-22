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
	
	public ObjectGoal<TaxTree, GSProject> getTaxTreeGoal() {
		return taxTreeGoal;
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
				getLogger().info("Number of completed tax ids: " + taxIdNodes.size());
			}
			
			set(taxIdNodes);
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}

}
