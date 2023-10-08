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