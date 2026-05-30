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
		if (getLogger().isDebugEnabled()) {
			logger.debug("Total added reads: " + added);
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