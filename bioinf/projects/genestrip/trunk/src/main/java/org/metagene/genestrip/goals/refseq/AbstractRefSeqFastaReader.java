package org.metagene.genestrip.goals.refseq;

import java.io.IOException;
import java.util.Set;

import org.metagene.genestrip.fasta.AbstractFastaReader;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;
import org.metagene.genestrip.util.ByteArrayUtil;

public abstract class AbstractRefSeqFastaReader extends AbstractFastaReader {
	private final Set<TaxIdNode> taxNodes;
	private final AccessionMap accessionMap;
	
	protected boolean includeRegion;
	protected TaxIdNode node;

	public AbstractRefSeqFastaReader(int bufferSize, Set<TaxIdNode> taxNodes, AccessionMap accessionMap) {
		super(bufferSize);
		this.taxNodes = taxNodes;
		this.accessionMap = accessionMap;
		includeRegion = false;
	}

	@Override
	protected void infoLine() throws IOException {
		includeRegion = false;
		int pos = ByteArrayUtil.indexOf(target, 0, size, ' ');
		if (pos >= 0) {
			node = accessionMap.get(target, 1, pos);
			if (node != null) {
				includeRegion = taxNodes.isEmpty() || taxNodes.contains(node);
			}
		}
	}
}
