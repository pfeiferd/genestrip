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
import org.metagene.genestrip.tax.SmallTaxTree;
import org.metagene.genestrip.tax.SmallTaxTree.SmallTaxIdNode;

import java.io.Serializable;
import java.util.Collections;
import java.util.Map;

public class MatchingResult implements Serializable {
	private static final long serialVersionUID = 1L;
	
	private final int k;
	private final Map<String, CountsPerTaxid> taxid2Stats;
	private final long totalReads;
	private final long totalKMers;
	private final short[] totalMaxCounts;

	public MatchingResult(int k, Map<String, CountsPerTaxid> taxid2Stats, long totalReads, long totalKMers,
			short[] totalMaxCounts) {
		this.k = k;
		this.taxid2Stats = taxid2Stats;
		this.totalReads = totalReads;
		this.totalKMers = totalKMers;
		this.totalMaxCounts = totalMaxCounts;
	}

	public void extendResults(Database database) {
		SmallTaxTree tree = database.getTaxTree();
		for (String key : taxid2Stats.keySet()) {
			SmallTaxIdNode node = tree.getNodeByTaxId(key);
			if (node == null) {
				for (node = node.getParent(); node != null; node = node.getParent()) {
					if (!taxid2Stats.containsKey(node.getTaxId())) {
						taxid2Stats.put(node.getTaxId(), new CountsPerTaxid(node.getTaxId(), 0));
					}
				}
			}
		}
		for (String key : taxid2Stats.keySet()) {
			CountsPerTaxid stats = taxid2Stats.get(key);
			if (stats != null) {
				long dbKMers = database.getStats().getLong(key);
				stats.setDbKMers(dbKMers);
				if (dbKMers > 0) {
					stats.setCoverage(((double) stats.getUniqueKMers()) / dbKMers);
					stats.setExpUnique(getExpectedUniqueKMers(stats));
					double nreads = ((double) stats.getReads()) / dbKMers;
					stats.setNormalizedReads(nreads);
					stats.addAccNormalizedReads(nreads);
					stats.addAccReads(stats.getReads());
					SmallTaxIdNode node = tree.getNodeByTaxId(key);
					if (node != null) {
						for (node = node.getParent(); node != null; node = node.getParent()) {
							stats = taxid2Stats.get(node.getTaxId());
							if (stats != null) {
								stats.addAccNormalizedReads(nreads);
								stats.addAccReads(stats.getReads());
							}
						}
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

	public long getTotalKMers() {
		return totalKMers;
	}

	public long getTotalReads() {
		return totalReads;
	}

	public short[] getTotalMaxCounts() {
		return totalMaxCounts;
	}

	public boolean isWithMaxKMerCounts() {
		return totalMaxCounts != null;
	}

	private double getExpectedUniqueKMers(CountsPerTaxid stats) {
		return  (1 - Math.pow(1 - 1d / stats.getDbKMers(), stats.getKMers())) * stats.getDbKMers();
	}
}