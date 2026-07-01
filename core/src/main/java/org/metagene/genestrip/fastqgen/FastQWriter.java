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
package org.metagene.genestrip.fastqgen;
import java.nio.charset.StandardCharsets;

import java.io.OutputStream;
import java.io.PrintStream;

import org.apache.commons.logging.Log;
import org.metagene.genestrip.util.GSLogFactory;

/**
 * Writes reads to an output stream in FASTQ format, generating a descriptor line (with an id and
 * optional tax id) and a placeholder quality line ('~') for each read.
 */
public class FastQWriter {
	private final Log logger = GSLogFactory.getLog("fastqwriter");

	private final String id;
	private PrintStream printStream;
	private long added;

	/**
	 * Creates a writer that emits reads with the given id prefix to the given output stream.
	 *
	 * @param id  the id prefix used in each read's descriptor line
	 * @param out the output stream to write the FASTQ data to
	 */
	public FastQWriter(String id, OutputStream out) {
		this.id = id;
		printStream = new PrintStream(out, false, StandardCharsets.UTF_8);
	}

	/**
	 * Returns the number of reads added so far.
	 *
	 * @return the number of reads added
	 */
	public long getAdded() {
		return added;
	}

	/**
	 * Returns the logger used by this writer.
	 *
	 * @return the logger
	 */
	protected Log getLogger() {
		return logger;
	}

	/**
	 * Writes one read consisting of the given base bytes, tagged with the given (optional) tax id.
	 *
	 * @param taxid the tax id to tag the read with, or {@code null} for none
	 * @param data  the read's base bytes
	 */
	public void addRead(String taxid, byte[] data) {
		added++;
		printBeforeRead(taxid);
		for (int i = 0; i < data.length; i++) {
			printStream.print((char) data[i]);
		}
		printStream.println();
		printAfterRead(data.length);
	}

	/**
	 * Finishes writing, logging the total number of reads added.
	 */
	public void done() {
		if (getLogger().isDebugEnabled()) {
			logger.debug("Total added reads: " + added);
		}
	}

	/**
	 * Resets the added-read counter to begin a new run.
	 */
	public void start() {
		added = 0;
	}

	/**
	 * Writes a read's descriptor line: the id, the optional tax id and the running read number.
	 *
	 * @param taxid the tax id to include, or {@code null} to omit it
	 */
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


	/**
	 * Writes a read's separator ('+') line and a placeholder quality line of '~' of the given length.
	 *
	 * @param readLength the length of the read (number of quality characters to write)
	 */
	protected void printAfterRead(int readLength) {
		printStream.println("+");
		for (int i = 0; i < readLength; i++) {
			printStream.print('~');
		}
		printStream.println();
	}
}