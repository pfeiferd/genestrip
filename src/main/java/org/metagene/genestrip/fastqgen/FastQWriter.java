package org.metagene.genestrip.fastqgen;

import java.io.OutputStream;
import java.io.PrintStream;

import org.apache.commons.logging.Log;
import org.metagene.genestrip.util.GSLogFactory;

public class FastQWriter {
	private final Log logger = GSLogFactory.getLog("fastqwriter");

	private final String id;
	private PrintStream printStream;
	private long added;

	public FastQWriter(String id, OutputStream out) {
		this.id = id;
		printStream = new PrintStream(out);
	}

	public long getAdded() {
		return added;
	}

	protected Log getLogger() {
		return logger;
	}

	public void addRead(String taxid, byte[] data) {
		added++;
		printBeforeRead(taxid);
		for (int i = 0; i < data.length; i++) {
			printStream.print((char) data[i]);
		}
		printStream.println();
		printAfterRead(data.length);
	}

	public void done() {
		if (getLogger().isInfoEnabled()) {
			logger.info("Total added reads: " + added);
		}
	}

	public void start() {
		added = 0;
	}

	protected void printBeforeRead(String taxid) {
		printStream.print(id);
		if (taxid != null) {
			printStream.print(':');
			printStream.print(taxid);
		}
		printStream.print(':');
		printStream.print(added);
		printStream.println();
	}


	protected void printAfterRead(int readLength) {
		printStream.println("+");
		for (int i = 0; i < readLength; i++) {
			printStream.print('~');
		}
		printStream.println();
	}
}