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
import org.metagene.genestrip.io.BufferedLineReader;
import org.metagene.genestrip.io.StreamProvider;
import org.metagene.genestrip.util.GSLogFactory;

public abstract class AbstractFastaReader {
	private final Log logger = GSLogFactory.getLog("fastareader");
	
	private final BufferedLineReader bufferedLineReader;
	protected final byte[] target;
	protected int size;
	
	protected long dataLines;
	
	public AbstractFastaReader(int bufferSize) {
		this.bufferedLineReader = new BufferedLineReader();
		target = new byte[bufferSize];
	}

	public void readFasta(File file) throws IOException {
		if (logger.isDebugEnabled()) {
			logger.debug("Reading fasta file " + file);
		}
		try (InputStream inputStream = StreamProvider.getInputStreamForFile(file)) {
			readFasta(inputStream);
		}
		if (logger.isDebugEnabled()) {
			logger.debug("Closed fasta file " + file);
		}
	}
	
	protected Log getLogger() {
		return logger;
	}

	// value must not be null...
	public void readFasta(InputStream inputStream) throws IOException {
		dataLines = 0;
		bufferedLineReader.setInputStream(inputStream);
		
		start();
		boolean first = true;
		for (size = bufferedLineReader.nextLine(target); size > 0; size = bufferedLineReader.nextLine(target)) {
			if (size >= target.length - 1) {
				throw new IllegalStateException("buffer is too small for data line in fastq file");
			}
			else {
				target[size] = 0;
			}
			if (target[0] == '>') {
				if (!first) {
					endRegion();
				}
				first = false;
				startRegion();
				infoLine();
			}
			else {
				dataLines++;
				dataLine();
			}
		}
		if (!first) {
			endRegion();
		}
		done();
		if (logger.isDebugEnabled()) {
			logger.debug("Total number of data lines: " + dataLines);
		}
	}

	protected void startRegion() {
	}

	protected void endRegion() {
	}
	
	protected void infoLine() {
	}
	
	protected abstract void dataLine();
	
	protected void done() {};
	
	protected void start() {};
}