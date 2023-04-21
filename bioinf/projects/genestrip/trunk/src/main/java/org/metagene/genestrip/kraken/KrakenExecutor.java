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
package org.metagene.genestrip.kraken;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.text.MessageFormat;

import org.apache.commons.io.IOUtils;
import org.metagene.genestrip.util.StreamProvider;

public class KrakenExecutor {
	private final String binFolder;
	private final String execCommand;

	public KrakenExecutor(String binFolder, String execCommand) {
		this.execCommand = execCommand;
		this.binFolder = binFolder;
	}

	public String genExecLine(String database, File fastq) throws IOException {
		return MessageFormat.format(execCommand, binFolder, database, fastq.getCanonicalPath());
	}

	public void execute(String database, File fastq, File outFile) throws InterruptedException, IOException {
		Process process = Runtime.getRuntime().exec(genExecLine(database, fastq));
		handleOutputStream(process.getInputStream(), outFile);
		int res = process.waitFor();
		if (res != 0) {
			throw new IllegalStateException("Kraken terminated unsuccesfully");
		}
	}

	protected void handleOutputStream(InputStream stream, File outFile) throws IOException {
		OutputStream out = StreamProvider.getOutputStreamForFile(outFile);
		IOUtils.copyLarge(stream, out);
		out.close();
	}
}
