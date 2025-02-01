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
package org.metagene.genestrip.refseq;

import java.util.Set;

import org.metagene.genestrip.fasta.AbstractFastaReader;
import org.metagene.genestrip.tax.Rank;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;
import org.metagene.genestrip.util.ByteArrayUtil;
import org.metagene.genestrip.util.StringLongDigitTrie;

public abstract class AbstractRefSeqFastaReader extends AbstractFastaReader {
	protected final Set<TaxIdNode> taxNodes;
	protected final AccessionMap accessionMap;
	protected final int maxGenomesPerTaxId;
	protected final long maxKmersPerTaxId;
	protected final Rank maxGenomesPerTaxIdRank;
	protected final StringLong2DigitTrie regionsPerTaxid;
	protected final int k;

	protected boolean includeRegion;
	protected long kMersForNode;
	protected TaxIdNode node;

	protected boolean ignoreMap;
	protected long includedCounter;
	private long bpsInRegion;
	protected long kmersInRegion;
	protected long totalKmers;

	public AbstractRefSeqFastaReader(int bufferSize, Set<TaxIdNode> taxNodes, AccessionMap accessionMap, int k, int maxGenomesPerTaxId, Rank maxGenomesPerTaxIdRank, long maxKmersPerTaxId) {
		super(bufferSize);
		this.taxNodes = taxNodes;
		this.accessionMap = accessionMap;
		this.k = k;
		includeRegion = false;
		ignoreMap = false;
		includedCounter = 0;
		totalKmers = 0;
		regionsPerTaxid = new StringLong2DigitTrie();
		this.maxGenomesPerTaxId = maxGenomesPerTaxId;
		this.maxGenomesPerTaxIdRank = maxGenomesPerTaxIdRank;
		this.maxKmersPerTaxId = maxKmersPerTaxId;
	}
	
	public StringLongDigitTrie getRegionsPerTaxid() {
		return regionsPerTaxid;
	}

	public void ignoreAccessionMap(TaxIdNode node) {
		this.ignoreMap = node != null;
		this.node = node;
	}

	@Override
	protected void startRegion() {
		includeRegion = false;
		kmersInRegion = 0;
		bpsInRegion = 0;
	}

	@Override
	protected void endRegion() {
		if (includeRegion) {
			totalKmers += kmersInRegion;
			if (node != null) {
				for (TaxIdNode n = node; n != null; n = n.getParent()) {
					regionsPerTaxid.incAndAdd(n.getTaxId(), kmersInRegion);
				}
			}
			includedCounter++;
		}
	}

	@Override
	protected void infoLine() {
		// Handle new region:
		if (!ignoreMap) {
			updateNodeFromInfoLine();
		}
		if (node != null && (taxNodes.isEmpty() || taxNodes.contains(node))) {
			includeRegion = true;
			kMersForNode = 0;
			if (maxGenomesPerTaxIdRank == null) {
				StringLong2DigitTrie.StringLong2 sl = (StringLong2DigitTrie.StringLong2) regionsPerTaxid.get(node.getTaxId());
				if (sl != null) {
					kMersForNode = sl.longValue2;
					if (kMersForNode >= maxKmersPerTaxId || sl.getLongValue() >= maxGenomesPerTaxId) {
						includeRegion = false;
					}
				}
			}
			else {
				for (TaxIdNode n = node; n != null; n = n.getParent()) {
					if (maxGenomesPerTaxIdRank.equals(n.getRank())) {
						StringLong2DigitTrie.StringLong2 sl =
								(StringLong2DigitTrie.StringLong2) regionsPerTaxid.get(n.getTaxId());
						if (sl != null) {
							kMersForNode = sl.longValue2;
							if (kMersForNode >= maxKmersPerTaxId || sl.getLongValue() >= maxGenomesPerTaxId) {
								includeRegion = false;
							}
						}
						break;
					}
				}
			}
		}
		else {
			includeRegion = false;
		}
	}

	@Override
	protected void dataLine() {
		if (includeRegion) {
			if (isAllowMoreKmers()) {
				bpsInRegion += size - 1;
				if (bpsInRegion >= k) {
					kmersInRegion = bpsInRegion - k + 1;
				}
			}
		}
	}

	public boolean isAllowMoreKmers() {
		return kMersForNode + kmersInRegion < maxKmersPerTaxId;
	}

	protected void updateNodeFromInfoLine() {
		int pos = ByteArrayUtil.indexOf(target, 0, size, ' ');
		if (pos >= 0) {
			node = accessionMap.get(target, 1, pos);
		}
		else {
			if (getLogger().isInfoEnabled()) {
				getLogger().info("Inconsistent info line: " + new String(target, 0, size - 1));
			}
			node = null;
		}
	}
	
	@Override
	protected void done() {
		super.done();
		if (getLogger().isInfoEnabled()) {
			getLogger().info("Number of included regions: " + includedCounter);
			getLogger().info("Total included kmers: " + totalKmers);
			getLogger().info("Resulting approx. DB Size in MB (without Bloom filter): " + (totalKmers * 10) / (1024 * 1024) );
		}
	}

	protected static class StringLong2DigitTrie extends StringLongDigitTrie {
		public void incAndAdd(String key, long add) {
			StringLong2 stringLong = (StringLong2) inc(key);
			stringLong.longValue2 += add;
		}

		@Override
		protected StringLong2 createInGet(String digits, Object createContext) {
			return new StringLong2(digits);
		}

		@Override
		protected StringLong createInGet(byte[] seq, int start, int end, Object createContext) {
			return new StringLong2(new String(seq, start, end - start));
		}


		public static class StringLong2 extends StringLong {
			private long longValue2;

			public StringLong2(String digits) {
				super(digits);
			}

			@Override
			public String toString() {
				return "SL:(taxid:" + stringValue + ", regions: " + longValue + ", kmers: " + longValue2 + ")";
			}
		}
	}
}
