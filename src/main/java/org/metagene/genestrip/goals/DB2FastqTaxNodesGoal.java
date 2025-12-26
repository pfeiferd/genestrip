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

import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import org.metagene.genestrip.GSConfigKey;
import org.metagene.genestrip.GSGoalKey;
import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.store.Database;
import org.metagene.genestrip.store.KMerSortedArray;
import org.metagene.genestrip.tax.Rank;
import org.metagene.genestrip.tax.SmallTaxTree;
import org.metagene.genestrip.tax.SmallTaxTree.SmallTaxIdNode;

public class DB2FastqTaxNodesGoal extends ObjectGoal<Set<SmallTaxIdNode>, GSProject> {
	private final ObjectGoal<Database, GSProject> storeGoal;

	@SafeVarargs
	public DB2FastqTaxNodesGoal(GSProject project,
			ObjectGoal<Database, GSProject> storeGoal, Goal<GSProject>... deps) {
		super(project, GSGoalKey.DB2FASTQ_TAXIDS, Goal.append(deps, storeGoal));
		this.storeGoal = storeGoal;
	}

	@Override
	protected void doMakeThis() {
		Set<SmallTaxIdNode> taxIdNodes;
		KMerSortedArray<String> store = storeGoal.get().getKmerStore();
		SmallTaxTree taxTree = storeGoal.get().getTaxTree();

		Set<SmallTaxIdNode> storedNodes = new HashSet<SmallTaxIdNode>();
		Iterator<String> it = store.getValues();
		while (it.hasNext()) {
			String storeId = it.next();
			if (storeId != null) {
				SmallTaxIdNode node = taxTree.getNodeByTaxId(storeId);
				if (node != null) {
					storedNodes.add(node);
				}
			}
		}

		List<String> taxidList = (List<String>)getProject().configValue(GSConfigKey.TAX_IDS);
		if (taxidList != null && !taxidList.isEmpty()) {
			String[] taxids = taxidList.toArray(new String[taxidList.size()]);
			boolean[] withDescs = new boolean[taxids.length];
			for (int i = 0; i < taxids.length; i++) {
				if (taxids[i].endsWith("+")) {
					taxids[i] = taxids[i].substring(0, taxids[i].length() - 1);
					withDescs[i] = true;
				}
			}

			taxIdNodes = asNodesWithDesc(taxTree, taxids, withDescs);
			taxIdNodes.retainAll(storedNodes);
		}
		else {
			taxIdNodes = storedNodes;
		}
		if (getLogger().isDebugEnabled()) {
			getLogger().debug("Requested tax ids: " + taxIdNodes);
		}
		set(taxIdNodes);
	}

	public Set<SmallTaxIdNode> asNodesWithDesc(SmallTaxTree taxTree, String[] taxids, boolean[] withDescs) {
		Set<SmallTaxIdNode> res = new HashSet<SmallTaxIdNode>();
		for (int i = 0; i < taxids.length; i++) {
			SmallTaxIdNode node = taxTree.getNodeByTaxId(taxids[i]);
			if (node != null) {
				res.add(node);
				if (withDescs[i]) {
					completeFilterlist(res, node, null);
				}
			}
		}

		return res;
	}

	private void completeFilterlist(Set<SmallTaxIdNode> filter, SmallTaxIdNode node, Rank depth) {
		if (node != null) {
			filter.add(node);
			SmallTaxIdNode[] subNodes = node.getSubNodes();
			if (subNodes != null) {
				for (SmallTaxIdNode subNode : subNodes) {
					if (depth == null
							|| (subNode != null && subNode.getRank() != null && !subNode.getRank().isBelow(depth))) {
						completeFilterlist(filter, subNode, depth);
					}
				}
			}
		}
	}	
}
