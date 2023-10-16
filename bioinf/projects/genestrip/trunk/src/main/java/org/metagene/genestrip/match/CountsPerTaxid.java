package org.metagene.genestrip.match;

import java.io.Serializable;

public class CountsPerTaxid implements Serializable {
	private static final long serialVersionUID = 1L;

	protected String taxid;
	protected long reads; // Only classified reads are counted - kraken style classification approach...
	protected long uniqueKmers; // All unique kmers counted - even from unclassified reads.
	protected long kmers; // All kmers counted - even from unclassified reads.
	protected int maxContigLen;
	protected int contigs;
	protected short[] maxKMerCounts;
	protected byte[] maxContigDescriptor;

	public CountsPerTaxid(String taxid, int maxReadSizeBytes) {
		this.taxid = taxid;
		maxContigDescriptor = new byte[maxReadSizeBytes];
	}

	public String getTaxid() {
		return taxid;
	}

	public int getContigs() {
		return contigs;
	}

	public long getKMers() {
		return kmers;
	}

	public int getMaxContigLen() {
		return maxContigLen;
	}

	public byte[] getMaxContigDescriptor() {
		return maxContigDescriptor;
	}

	public long getReads() {
		return reads;
	}

	public long getUniqueKMers() {
		return uniqueKmers;
	}

	public short[] getMaxKMerCounts() {
		return maxKMerCounts;
	}
}