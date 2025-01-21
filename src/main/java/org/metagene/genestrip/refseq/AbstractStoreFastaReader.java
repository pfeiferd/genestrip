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

import org.metagene.genestrip.tax.Rank;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;
import org.metagene.genestrip.util.CGAT;
import org.metagene.genestrip.util.CGATLongBuffer;
import org.metagene.genestrip.util.CGATRingBuffer;

public abstract class AbstractStoreFastaReader extends AbstractRefSeqFastaReader {
	protected final CGATLongBuffer byteRingBuffer;
	protected long dustCounter;
	protected long totalKmerCounter;
	
	public AbstractStoreFastaReader(int bufferSize, Set<TaxIdNode> taxNodes, AccessionMap accessionMap, int k, int maxGenomesPerTaxId, Rank maxGenomesPerTaxIdRank, long maxKmersPerTaxId, int maxDust) {
		super(bufferSize, taxNodes, accessionMap, maxGenomesPerTaxId, maxGenomesPerTaxIdRank, maxKmersPerTaxId);
		byteRingBuffer = k > 32 ? new CGATRingBuffer(k, maxDust) : new CGATLongBuffer(k, maxDust);
		dustCounter = 0;
	}
	
	@Override
	protected void infoLine() throws IOException {
		super.infoLine();
		byteRingBuffer.reset();
	}
	
	@Override
	protected void dataLine() throws IOException {
		if (includeRegion) {
			for (int i = 0; i < size - 1; i++) {
				byteRingBuffer.put(CGAT.cgatToUpperCase(target[i]));
				if (byteRingBuffer.isFilled()) {
					 if (!byteRingBuffer.isDust()) {
						 handleStore();
						 kmersInRegion++;
					 }
					 else {
						 dustCounter++;
					 }
					 totalKmerCounter++;
				}
			}
		}
	}
	
	@Override
	protected void done() throws IOException {
		super.done();
		if (getLogger().isInfoEnabled()) {
			getLogger().info("Dust ratio: " + ((double) dustCounter) / totalKmerCounter);
		}
	}
	
	protected abstract void handleStore();
}
