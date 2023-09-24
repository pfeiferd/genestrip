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
package org.metagene.genestrip.io;

import java.io.IOException;
import java.io.InputStream;

public class BufferedLineReader {
	public static final int DEFAULT_BUFFER_SIZE = 8 * 4096;
	
	private static int bufferSize = DEFAULT_BUFFER_SIZE;
	
	public static void setBufferSize(int bufferSize) {
		BufferedLineReader.bufferSize = bufferSize;
	}
	
	public static int getBufferSize() {
		return bufferSize;
	}
	
	private final byte[] buffer;

	private int pos; // Position from which to start reading from the buffer.
	private InputStream stream;
	private int bufferFill;
	
	private long millis;

	public BufferedLineReader() {
		this(null);
	}
	
	public BufferedLineReader(InputStream stream) {
		this(stream, bufferSize);
	}
	
	public BufferedLineReader(InputStream stream, int bufferSize) {
		this.buffer = new byte[bufferSize];
		setInputStream(stream);
	}

	public void setInputStream(InputStream stream) {
		pos = buffer.length;
		bufferFill = 0;
		this.stream = stream;
	}

	public int skipLine() throws IOException {
		int size = 0;
		while (bufferFill != -1) {
			if (pos < bufferFill) {
				byte c = -1;
				for (; pos < bufferFill && c != '\n'; pos++) {
					c = buffer[pos];
					if (c != 0) {
						size++;
					}
				}
				if (c == '\n') {
					return size;
				}
			}
			long start = System.currentTimeMillis();
			bufferFill = stream.read(buffer);
			millis += System.currentTimeMillis() - start;
			pos = 0;
		}
		return size;
	}
	
	public int nextLine(byte[] target) throws IOException {
		int size = 0;
		while (bufferFill != -1) {
			if (pos < bufferFill) {
				byte c = -1;
				for (; size < target.length && pos < bufferFill && c != '\n'; pos++) {
					target[size] = c = buffer[pos];
					if (c != 0) {
						size++;
					}
				}
				if (size == target.length || c == '\n') {
					return size;
				}
			}
			long start = System.currentTimeMillis();
			bufferFill = stream.read(buffer);
			millis += System.currentTimeMillis() - start;
			pos = 0;
		}
		return size;
	}
	
	public long getMillis() {
		return millis;
	}
	
	public void close() throws IOException {
		stream.close();
	}
}
