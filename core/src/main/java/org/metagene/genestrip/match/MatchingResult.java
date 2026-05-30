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
package org.metagene.genestrip.match;

import org.metagene.genestrip.store.Database;
import org.metagene.genestrip.tax.Rank;
import org.metagene.genestrip.tax.SmallTaxTree;
import org.metagene.genestrip.tax.SmallTaxTree.SmallTaxIdNode;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;

public class MatchingResult implements Serializable {
	private static final long serialVersionUID = 1L;
	
	private final int k;
	private final Map<String, CountsPerTaxid> taxid2Stats;
	private final CountsPerTaxid globalStats;

	public MatchingResult(int k, Map<String, CountsPerTaxid> taxid2Stats, long totalReads, long totalKMers, long totalBPs,
			short[] totalMaxCounts) {
		this.k = k;
		globalStats =  new CountsPerTaxid(0, null, totalReads, totalKMers, totalBPs, totalMaxCounts);
		this.taxid2Stats = taxid2Stats;
	}


	public void completeResults(Database database) {
		taxid2Stats.put(globalStats.getTaxid(), globalStats);

		SmallTaxTree tree = database.getTaxTree();
		// Add potentially missing parent nodes for report.
		for (String key : new ArrayList<>(taxid2Stats.keySet())) {
			SmallTaxIdNode node = tree.getNodeByTaxId(key);
			if (node != null) {
				for (node = node.getParent(); node != null; node = node.getParent()) {
					if (!taxid2Stats.containsKey(node.getTaxId())) {
						taxid2Stats.put(node.getTaxId(), new CountsPerTaxid(node.getLevel(), node.getTaxId(), 0));
					}
				}
			}
		}

		//globalStats.completeValues(pos++, database.getStats().getLong(null), tree.getNodeByTaxId("1"));
		List<String> keys = new ArrayList<String>(taxid2Stats.keySet());
		tree.sortTaxidsViaTree(keys);
		int pos = 0;
		for (String key : keys) {
			CountsPerTaxid stats = taxid2Stats.get(key);
			long dbKMers = database.getStats().getLong(key);
			SmallTaxIdNode node = tree.getNodeByTaxId(key);
			stats.completeValues(pos++, dbKMers, node);
			if (node != null) {
				for (node = node.getParent(); node != null; node = node.getParent()) {
					CountsPerTaxid stats2 = taxid2Stats.get(node.getTaxId());
					if (stats2 != null) {
						stats2.accumulateFrom(stats);
					}
				}
			}
		}
	}

	public int getK() {
		return k;
	}

	public Map<String, CountsPerTaxid> getTaxid2Stats() {
		return Collections.unmodifiableMap(taxid2Stats);
	}

	public boolean isWithMaxKMerCounts() {
		return globalStats.maxKMerCounts != null;
	}

	public CountsPerTaxid getGlobalStats() {
		return globalStats;
	}
}