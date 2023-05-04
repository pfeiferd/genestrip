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

import org.apache.commons.exec.CommandLine;
import org.apache.commons.exec.util.StringUtils;
import org.apache.commons.io.IOUtils;

public class SortExecutor {
	private final String bin;
	private final String execCommand;

	public SortExecutor(String bin, String execCommand) {
		this.execCommand = execCommand;
		this.bin = StringUtils.quoteArgument(bin);
	}

	public String genExecLine(File fastqIn, File fastqOut) throws IOException {
		return MessageFormat.format(execCommand, bin, StringUtils.quoteArgument(fastqIn.getCanonicalPath()),
				StringUtils.quoteArgument(fastqOut.getCanonicalPath()));
	}

	protected CommandLine genExecCommand(File fastqIn, File fastqOut) throws IOException {
		return CommandLine.parse(genExecLine(fastqIn, fastqOut));
	}

	public void execute(File fastqIn, File fastqOut)
			throws InterruptedException, IOException {
		Process process = Runtime.getRuntime().exec(genExecLine(fastqIn, fastqOut));
		handleErrorStream(process.getErrorStream(), System.err);
		int res = process.waitFor();
		if (res != 0) {
			throw new IllegalStateException("Sort terminated unsuccesfully");
		}
	}

	protected void handleErrorStream(InputStream inputStream, OutputStream errorStream) throws IOException {
		IOUtils.copy(inputStream, errorStream);
	}
}
