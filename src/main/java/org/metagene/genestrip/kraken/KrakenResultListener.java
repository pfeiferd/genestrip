package org.metagene.genestrip.kraken;

public interface KrakenResultListener {
	public void newClassInfo(long lineCount, String taxid, int fr, byte[] readDescriptor);
}
