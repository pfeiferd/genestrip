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

import org.metagene.genestrip.tax.Rank;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;
import org.metagene.genestrip.util.CGAT;
import org.metagene.genestrip.util.CGATLongBuffer;

public abstract class AbstractStoreFastaReader extends AbstractRefSeqFastaReader {
	protected final CGATLongBuffer byteRingBuffer;
	protected long dustCounter;
	protected long totalKmers;

	private boolean enableLowerCaseBases;

	public AbstractStoreFastaReader(int bufferSize, Set<TaxIdNode> taxNodes, AccessionMap accessionMap, int k, int maxGenomesPerTaxId, Rank maxGenomesPerTaxIdRank,
									long maxKmersPerTaxId, int maxDust, int stepSize, boolean completeGenomesOnly, StringLong2DigitTrie regionsPerTaxid,
									boolean enableLowerCaseBases) {
		super(bufferSize, taxNodes, accessionMap, k, maxGenomesPerTaxId, maxGenomesPerTaxIdRank, maxKmersPerTaxId, stepSize, completeGenomesOnly, regionsPerTaxid);
		byteRingBuffer = new CGATLongBuffer(k, maxDust);
		dustCounter = 0;
		this.enableLowerCaseBases = enableLowerCaseBases;
	}

	@Override
	protected void startRegion() {
		super.startRegion();
		byteRingBuffer.reset();
	}

	@Override
	protected void dataLine() {
		if (includeRegion) {
			if (isAllowMoreKmers()) {
				for (int i = 0; i < size - 1; i++) {
					byteRingBuffer.put(enableLowerCaseBases ? CGAT.cgatToUpperCase(target[i]) : target[i]);
					bpsInRegion++;
					if (bpsInRegion % stepSize == 0) {
						if (byteRingBuffer.isFilled()) {
							if (byteRingBuffer.isDust()) {
								dustCounter++;
							} else if (handleStore()) {
								kmersInRegion++;
							}
							totalKmers++;
						}
					}
				}
			}
		}
	}

	protected abstract boolean handleStore();
}
