package org.metagene.genestrip.goals.refseq;

import org.metagene.genestrip.tax.TaxTree.TaxIdNode;

public interface AccessionMap {
	public void put(byte[] array, int start, int end, TaxIdNode node);
	public TaxIdNode get(byte[] array, int start, int end);
	public void optimize();
}
