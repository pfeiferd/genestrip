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

public class CountsPerTaxid implements Serializable {
	private static final long serialVersionUID = 1L;

	protected String taxid;
	protected long reads; // Only classified if a read is ENTIRELY consistent with the kmers of a taxon (and/or its ancestors) in the DB
	protected long readKmers; // Kmers from classified reads.
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
	
	public long getReadKMers() {
		return readKmers;
	}

	public long getUniqueKMers() {
		return uniqueKmers;
	}

	public short[] getMaxKMerCounts() {
		return maxKMerCounts;
	}
}