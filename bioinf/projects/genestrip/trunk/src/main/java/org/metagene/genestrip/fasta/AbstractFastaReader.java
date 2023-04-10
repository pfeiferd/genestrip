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
package org.metagene.genestrip.fasta;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.metagene.genestrip.util.BufferedLineReader;
import org.metagene.genestrip.util.StreamProvider;

public abstract class AbstractFastaReader {
	private final Log logger = LogFactory.getLog(getClass());
	
	private final BufferedLineReader bufferedLineReader;
	protected final byte[] target;
	protected int size;
	
	public AbstractFastaReader(int bufferSize) {
		this.bufferedLineReader = new BufferedLineReader();
		target = new byte[bufferSize];
	}

	public void readFasta(File file) throws IOException {
		if (logger.isInfoEnabled()) {
			logger.info("Reading fasta file " + file);
		}
		InputStream inputStream = StreamProvider.getInputStreamForFile(file);
		readFasta(inputStream);
		inputStream.close();
		if (logger.isInfoEnabled()) {
			logger.info("Closed fasta file " + file);
		}
	}
	
	public Log getLogger() {
		return logger;
	}

	// value must not be null...
	public void readFasta(InputStream inputStream) throws IOException {
		bufferedLineReader.setInputStream(inputStream);
		
		for (size = bufferedLineReader.nextLine(target); size > 0; size = bufferedLineReader.nextLine(target)) {
			if (size >= target.length - 1) {
				throw new IllegalStateException("buffer is too small for data line in fastq file");
			}
			else {
				target[size] = 0;
			}
			if (target[0] == '>') {
				infoLine();
			}
			else {
				dataLine();
			}
		}
	}
	
	protected void infoLine() {
	}
	
	protected abstract void dataLine();
}