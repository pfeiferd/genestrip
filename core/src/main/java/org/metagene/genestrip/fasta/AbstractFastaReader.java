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

/**
 * Abstract line-based FASTA reader that reads a stream line by line and dispatches header ('&gt;')
 * lines and data (sequence) lines to overridable template methods.
 */
public abstract class AbstractFastaReader {
	private final Log logger = GSLogFactory.getLog("fastareader");
	
	private final BufferedLineReader bufferedLineReader;
	/** Buffer holding the bytes of the line currently being processed. */
	protected final byte[] target;
	/** Number of valid bytes of the current line held in {@link #target}. */
	protected int size;

	/** Number of data (sequence) lines read so far from the current input. */
	protected long dataLines;
	
	/**
	 * Creates a reader using a line buffer of the given size in bytes.
	 *
	 * @param bufferSize the size in bytes of the line buffer used to hold each line.
	 */
	public AbstractFastaReader(int bufferSize) {
		this.bufferedLineReader = new BufferedLineReader();
		target = new byte[bufferSize];
	}

	/**
	 * Reads the given FASTA file (transparently decompressing it as needed).
	 *
	 * @param file the FASTA file to read.
	 * @throws IOException if the file cannot be opened or read.
	 */
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
	
	/**
	 * Returns the logger used by this reader.
	 *
	 * @return the logger of this reader.
	 */
	protected Log getLogger() {
		return logger;
	}

	/**
	 * Reads FASTA content from the given (non-{@code null}) stream, dispatching header and data lines
	 * to the overridable template methods.
	 *
	 * @param inputStream the stream to read FASTA content from; must not be {@code null}.
	 * @throws IOException if the stream cannot be read.
	 */
	// value must not be null...
	public void readFasta(InputStream inputStream) throws IOException {
		dataLines = 0;
		bufferedLineReader.setInputStream(inputStream);
		
		start();
		boolean first = true;
		for (size = bufferedLineReader.nextLine(target); size > 0; size = bufferedLineReader.nextLine(target)) {
			if (size >= target.length - 1) {
				throw new IllegalStateException("buffer is too small for data line in fasta file");
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

	/**
	 * Called when a new region begins (i.e. a '&gt;' header line is encountered).
	 */
	protected void startRegion() {
	}

	/**
	 * Called when the current region ends (before the next region or at end of input).
	 */
	protected void endRegion() {
	}
	
	/**
	 * Called for a header ('&gt;') line; the line's bytes are held in {@code target[0, size)}.
	 */
	protected void infoLine() {
	}
	
	/**
	 * Handles a data (sequence) line, whose bytes are held in {@code target[0, size)}.
	 */
	protected abstract void dataLine();
	
	/**
	 * Called once after all lines of the input have been processed.
	 */
	protected void done() {};

	/**
	 * Called once before reading of the input begins.
	 */
	protected void start() {};
}