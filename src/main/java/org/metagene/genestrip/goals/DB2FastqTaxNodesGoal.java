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
import java.util.Set;

import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.store.KMerSortedArray;
import org.metagene.genestrip.store.KMerStoreWrapper;
import org.metagene.genestrip.tax.TaxIdCollector;
import org.metagene.genestrip.tax.TaxTree;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;

public class DB2FastqTaxNodesGoal extends ObjectGoal<Set<TaxIdNode>, GSProject> {
	private final ObjectGoal<TaxTree, GSProject> taxTreeGoal;
	private final ObjectGoal<KMerStoreWrapper, GSProject> storeGoal;

	@SafeVarargs
	public DB2FastqTaxNodesGoal(GSProject project, String name, ObjectGoal<TaxTree, GSProject> taxTreeGoal,
			ObjectGoal<KMerStoreWrapper, GSProject> storeGoal, Goal<GSProject>... deps) {
		super(project, name, Goal.append(deps, taxTreeGoal, storeGoal));
		this.taxTreeGoal = taxTreeGoal;
		this.storeGoal = storeGoal;
	}

	@Override
	public void makeThis() {
		Set<TaxIdNode> taxIdNodes;
		TaxTree taxTree = taxTreeGoal.get();
		KMerSortedArray<String> store = storeGoal.get().getKmerStore();

		Set<TaxIdNode> storedNodes = new HashSet<TaxTree.TaxIdNode>();
		for (short j = 0; j < store.getNValues(); j++) {
			String storeId = store.getValueForIndex(j);
			if (storeId != null) {
				TaxIdNode node = taxTree.getNodeByTaxId(storeId);
				if (node != null) {
					storedNodes.add(node);						
				}
			}
		}

		String encodedTaxIds = getProject().getTaxids();
		if (encodedTaxIds != null) {
			TaxIdCollector taxIdCollector = new TaxIdCollector(taxTree);
			String[] taxids = encodedTaxIds.split(",");
			boolean[] withDescs = new boolean[taxids.length];
			for (int i = 0; i < taxids.length; i++) {
				if (taxids[i].endsWith("+")) {
					taxids[i] = taxids[i].substring(0, taxids[i].length() - 1);
					withDescs[i] = true;
				}
			}

			taxIdNodes = taxIdCollector.asNodesWithDesc(taxids, withDescs);
			for (short j = 0; j < store.getNValues(); j++) {
				String storeId = store.getValueForIndex(j);
				if (storeId != null) {
				}
			}
			taxIdNodes.retainAll(storedNodes);
		}
		else {
			taxIdNodes = storedNodes;
		}
		if (getLogger().isInfoEnabled()) {
			getLogger().info("Requested tax ids: " + taxIdNodes);
		}
		set(taxIdNodes);
	}
}
