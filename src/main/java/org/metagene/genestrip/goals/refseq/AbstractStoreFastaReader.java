package org.metagene.genestrip.goals.refseq;

import java.io.IOException;
import java.util.Set;

import org.metagene.genestrip.tax.TaxTree.TaxIdNode;
import org.metagene.genestrip.util.CGAT;
import org.metagene.genestrip.util.CGATRingBuffer;

public abstract class AbstractStoreFastaReader extends AbstractRefSeqFastaReader {
	protected final CGATRingBuffer byteRingBuffer;
	
	public AbstractStoreFastaReader(int bufferSize, Set<TaxIdNode> taxNodes, AccessionMap accessionMap, int k) {
		super(bufferSize, taxNodes, accessionMap);
		byteRingBuffer = new CGATRingBuffer(k);
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
				if (byteRingBuffer.isFilled() && byteRingBuffer.isCGAT()) {
					handleStore(); 
				}
			}
		}
	}
	
	protected abstract void handleStore();
}
