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
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;
import org.metagene.genestrip.util.ByteArrayUtil;

public abstract class AbstractRefSeqFastaReader extends AbstractFastaReader {
	protected final Set<TaxIdNode> taxNodes;
	protected final AccessionMap accessionMap;

	protected boolean includeRegion;
	protected TaxIdNode node;

	protected boolean ignoreMap;

	public AbstractRefSeqFastaReader(int bufferSize, Set<TaxIdNode> taxNodes, AccessionMap accessionMap) {
		super(bufferSize);
		this.taxNodes = taxNodes;
		this.accessionMap = accessionMap;
		includeRegion = false;
	}

	public void ignoreAccessionMap(TaxIdNode node) {
		this.ignoreMap = true;
		this.node = node;
	}

	@Override
	protected void infoLine() throws IOException {
		if (!ignoreMap) {
			updateNodeFromInfoLine();
		}
		includeRegion = taxNodes.isEmpty() || (node != null && taxNodes.contains(node));
	}

	protected void updateNodeFromInfoLine() {
		int pos = ByteArrayUtil.indexOf(target, 0, size, ' ');
		if (pos >= 0) {
			node = accessionMap.get(target, 1, pos);
		}
	}
}
