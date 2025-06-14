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
	protected final int stepSize;
	private final boolean completeGenomesOnly;

	protected boolean includeRegion;
	protected long kMersForNode;
	protected TaxIdNode node;

	protected boolean ignoreMap;
	protected long bpsInRegion;
	protected long kmersInRegion;
	protected long includedKmers;


	public AbstractRefSeqFastaReader(int bufferSize, Set<TaxIdNode> taxNodes, AccessionMap accessionMap, int k, int maxGenomesPerTaxId,
									 Rank maxGenomesPerTaxIdRank, long maxKmersPerTaxId,
									 int stepSize, boolean completeGenomesOnly, StringLong2DigitTrie regionsPerTaxid) {
		super(bufferSize);
		this.taxNodes = taxNodes;
		this.accessionMap = accessionMap;
		this.k = k;
		this.stepSize = stepSize;
		includeRegion = false;
		ignoreMap = false;
		includedKmers = 0;
		this.regionsPerTaxid = regionsPerTaxid;
		this.maxGenomesPerTaxId = maxGenomesPerTaxId;
		this.maxGenomesPerTaxIdRank = maxGenomesPerTaxIdRank;
		this.maxKmersPerTaxId = maxKmersPerTaxId;
		this.completeGenomesOnly = completeGenomesOnly;
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
			includedKmers += kmersInRegion;
			if (node != null) {
				for (TaxIdNode n = node; n != null; n = n.getParent()) {
					regionsPerTaxid.incAndAdd(n.getTaxId(), kmersInRegion);
				}
			}
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

	public boolean isAllowMoreKmers() {
		return kMersForNode + kmersInRegion < maxKmersPerTaxId;
	}

	protected void updateNodeFromInfoLine() {
		int pos = ByteArrayUtil.indexOf(target, 0, size, ' ');
		if (pos >= 0) {
			node = accessionMap.get(target, 1, pos, completeGenomesOnly);
		}
		else {
			if (getLogger().isWarnEnabled()) {
				getLogger().warn("Inconsistent info line: " + new String(target, 0, size - 1));
			}
			node = null;
		}
	}

	public static class StringLong2DigitTrie extends StringLongDigitTrie {
		public void incAndAdd(String key, long add) {
			((StringLong2) get(key, this)).incAndAdd(add);
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

			public synchronized void incAndAdd(long add) {
				longValue++;
				longValue2 += add;
			}

			@Override
			public String toString() {
				return "SL:(taxid:" + stringValue + ", regions: " + longValue + ", kmers: " + longValue2 + ")";
			}
		}
	}
}
