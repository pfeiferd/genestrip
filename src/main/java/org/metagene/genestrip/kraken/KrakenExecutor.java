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
import org.apache.commons.exec.DefaultExecutor;
import org.apache.commons.exec.Executor;
import org.apache.commons.exec.PumpStreamHandler;
import org.apache.commons.exec.util.StringUtils;
import org.apache.commons.io.IOUtils;

public class KrakenExecutor {
	private final String bin;
	private final String execCommand;

	public KrakenExecutor(String bin, String execCommand) {
		this.execCommand = execCommand;
		this.bin = StringUtils.quoteArgument(bin);
	}

	public String genExecLine(String database, File fastq) throws IOException {
		return MessageFormat.format(execCommand, bin, StringUtils.quoteArgument(database),
				StringUtils.quoteArgument(fastq.getCanonicalPath()));
	}

	protected CommandLine genExecCommand(String database, File fastq) throws IOException {
		return CommandLine.parse(genExecLine(database, fastq));
	}

	public void execute2(String database, File fastq, OutputStream outputStream, OutputStream errorStream)
			throws IOException {
		Executor executor = new DefaultExecutor();
		PumpStreamHandler pumpStreamHandler = new PumpStreamHandler(outputStream, errorStream) {
			@Override
			public void setProcessOutputStream(final InputStream is) {
				// Event set it when it is null:
				createProcessOutputPump(is, outputStream);
			}

			protected Thread createPump(final InputStream is, final OutputStream os, final boolean closeWhenExhausted) {
				final Thread result = new Thread(new KrakenStreamPumper(is, os, closeWhenExhausted) {
					@Override
					protected void doWork(InputStream is, OutputStream os) throws IOException {
						if (os == outputStream) {
							handleOutputStream(is, os);
						} else if (os == errorStream) {
							handleErrorStream(is, os);
						} else {
							throw new IllegalStateException("bad stream in subsystem");
						}
					}
				}, "Kraken Exec Stream Pumper") {
				};
				result.setDaemon(true);
				return result;
			}
		};
		executor.setStreamHandler(pumpStreamHandler);
		int res = executor.execute(genExecCommand(database, fastq));
		if (res != 0) {
			throw new IllegalStateException("Kraken terminated unsuccesfully");
		}
	}

	public void execute(String database, File fastq, OutputStream outputStream, OutputStream errorStream)
			throws InterruptedException, IOException {
		Process process = Runtime.getRuntime().exec(genExecLine(database, fastq));
		handleOutputStream(process.getInputStream(), outputStream);
		handleErrorStream(process.getErrorStream(), System.err);
		int res = process.waitFor();
		if (res != 0) {
			throw new IllegalStateException("Kraken terminated unsuccesfully");
		}
	}

	protected void handleOutputStream(InputStream stream, OutputStream outputStream) throws IOException {
		IOUtils.copyLarge(stream, outputStream);
	}

	protected void handleErrorStream(InputStream inputStream, OutputStream errorStream) throws IOException {
		IOUtils.copy(inputStream, errorStream);
	}
}