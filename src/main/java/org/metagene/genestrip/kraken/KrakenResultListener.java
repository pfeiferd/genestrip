package org.metagene.genestrip.kraken;

public interface KrakenResultListener {
	public void newTaxIdForRead(long lineCount, byte[] readDescriptor, String krakenTaxid, int bps, String kmerTaxid, int hitLength);
}
