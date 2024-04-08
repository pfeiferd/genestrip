package org.metagene.genestrip.fastq;

public class FastqReaderInterruptedException extends RuntimeException {
	private static final long serialVersionUID = 1L;

	public FastqReaderInterruptedException(InterruptedException arg0) {
		super(arg0);
	}
}
