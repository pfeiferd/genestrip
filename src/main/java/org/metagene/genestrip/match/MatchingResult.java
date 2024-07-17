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
}