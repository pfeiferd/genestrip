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

import java.io.IOException;
import java.util.Set;

import org.metagene.genestrip.fasta.AbstractFastaReader;
import org.metagene.genestrip.tax.Rank;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;
import org.metagene.genestrip.util.ByteArrayUtil;
import org.metagene.genestrip.util.StringLongDigitTrie;
import org.metagene.genestrip.util.StringLongDigitTrie.StringLong;

public abstract class AbstractRefSeqFastaReader extends AbstractFastaReader {
	protected final Set<TaxIdNode> taxNodes;
	protected final AccessionMap accessionMap;
	protected final int maxGenomesPerTaxId;
	protected final Rank maxGenomesPerTaxIdRank;
	protected final StringLongDigitTrie regionsPerTaxid;

	protected boolean includeRegion;
	protected TaxIdNode node;

	protected boolean ignoreMap;
	protected long includedCounter;

	public AbstractRefSeqFastaReader(int bufferSize, Set<TaxIdNode> taxNodes, AccessionMap accessionMap, int maxGenomesPerTaxId, Rank maxGenomesPerTaxIdRank) {
		super(bufferSize);
		this.taxNodes = taxNodes;
		this.accessionMap = accessionMap;
		includeRegion = false;
		ignoreMap = false;
		includedCounter = 0;
		regionsPerTaxid = new StringLongDigitTrie();
		this.maxGenomesPerTaxId = maxGenomesPerTaxId;
		this.maxGenomesPerTaxIdRank = maxGenomesPerTaxIdRank;
	}
	
	public StringLongDigitTrie getRegionsPerTaxid() {
		return regionsPerTaxid;
	}

	public void ignoreAccessionMap(TaxIdNode node) {
		this.ignoreMap = node != null;
		this.node = node;
	}
	
	@Override
	protected void start() throws IOException {
		includeRegion = false;
	}

	@Override
	protected void infoLine() throws IOException {
		if (!ignoreMap) {
			updateNodeFromInfoLine();
		}
		if (node != null && (taxNodes.isEmpty() || taxNodes.contains(node))) {
			includeRegion = true;
			if (maxGenomesPerTaxIdRank == null) {
				StringLong sl = regionsPerTaxid.get(node.getTaxId());
				if (sl != null && sl.getLongValue() >= maxGenomesPerTaxId) {
					includeRegion = false;
				}
			}
			else {
				for (TaxIdNode n = node; n != null; n = n.getParent()) {
					if (maxGenomesPerTaxIdRank.equals(n.getRank())) {
						StringLong sl = regionsPerTaxid.get(n.getTaxId());
						if (sl != null && sl.getLongValue() >= maxGenomesPerTaxId) {
							includeRegion = false;
						}
						break;
					}
				}
			}
			if (includeRegion) {
				for (TaxIdNode n = node; n != null; n = n.getParent()) {
					regionsPerTaxid.inc(n.getTaxId());
				}
				includedCounter++;
			}
		}
		else {
			includeRegion = false;
		}
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
	protected void done() throws IOException {
		super.done();
		if (getLogger().isInfoEnabled()) {
			getLogger().info("Number of included regions: " + includedCounter);
		}
	}
}
