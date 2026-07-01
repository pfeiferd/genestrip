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

import java.io.Closeable;
import java.io.IOException;
import java.io.InputStream;

/**
 * Reads newline-delimited lines from an input stream into caller-supplied byte buffers with minimal
 * allocation. Null (0) bytes are skipped and not counted, and lines are delimited by {@code '\n'}.
 */
public class BufferedLineReader implements Closeable {
	/**
	 * Default internal buffer size (in bytes) used when no explicit size is given.
	 */
	public static final int DEFAULT_BUFFER_SIZE = 8 * 4096;

	private static int bufferSize = DEFAULT_BUFFER_SIZE;

	/**
	 * Sets the default internal buffer size (in bytes) used by newly created readers.
	 *
	 * @param bufferSize the default buffer size in bytes
	 */
	public static void setBufferSize(int bufferSize) {
		BufferedLineReader.bufferSize = bufferSize;
	}

	/**
	 * Returns the default internal buffer size (in bytes) used by newly created readers.
	 *
	 * @return the default buffer size in bytes
	 */
	public static int getBufferSize() {
		return bufferSize;
	}

	private final byte[] buffer;

	private int pos; // Position from which to start reading from the buffer.
	private InputStream stream;
	private int bufferFill;

	/**
	 * Creates a reader with no input stream set yet; use {@link #setInputStream(InputStream)} before
	 * reading.
	 */
	public BufferedLineReader() {
		this(null);
	}

	/**
	 * Creates a reader over the given stream using the current default buffer size.
	 *
	 * @param stream the input stream to read from
	 */
	public BufferedLineReader(InputStream stream) {
		this(stream, bufferSize);
	}

	/**
	 * Creates a reader over the given stream using an internal buffer of the given size.
	 *
	 * @param stream     the input stream to read from
	 * @param bufferSize the internal buffer size in bytes
	 */
	public BufferedLineReader(InputStream stream, int bufferSize) {
		this.buffer = new byte[bufferSize];
		setInputStream(stream);
	}

	/**
	 * Sets (or replaces) the input stream to read from and resets the internal buffer state.
	 *
	 * @param stream the input stream to read from
	 */
	public void setInputStream(InputStream stream) {
		pos = buffer.length;
		bufferFill = 0;
		this.stream = stream;
	}

	/**
	 * Reads and discards the next line, returning the number of non-null bytes it contained
	 * (excluding the terminating newline).
	 *
	 * @return the number of non-null bytes in the skipped line
	 * @throws IOException if reading from the underlying stream fails
	 */
	// Made final for potential inlining by JVM
	public final int skipLine() throws IOException {
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
			bufferFill = stream.read(buffer);
			pos = 0;
		}
		return size;
	}

	/**
	 * Reads the next line into {@code target} starting at position 0; see
	 * {@link #nextLine(byte[], int)}.
	 *
	 * @param target the buffer to write the line into
	 * @return the number of bytes written; {@code target.length + 1} if the buffer filled up before
	 *         the end of the line was reached
	 * @throws IOException if reading from the underlying stream fails
	 */
	// Made final for potential inlining by JVM
	public final int nextLine(byte[] target) throws IOException {
		return nextLine(target, 0);
	}

	/**
	 * Reads the next line into {@code target} starting at {@code startPos}, skipping null bytes.
	 *
	 * @param target   the buffer to write the line into
	 * @param startPos the position in {@code target} at which to start writing
	 * @return the number of bytes written; {@code target.length + 1} if the buffer filled up before
	 *         the end of the line was reached
	 * @throws IOException if reading from the underlying stream fails
	 */
	// Returns target.length + 1 if target is full but end of line not reached.
	// Made final for potential inlining by JVM
	public final int nextLine(byte[] target, int startPos) throws IOException {
		int size = startPos;
		while (bufferFill != -1) {
			if (pos < bufferFill) {
				byte c = -1;
				for (; size < target.length && pos < bufferFill && c != '\n'; pos++) {
					target[size] = c = buffer[pos];
					if (c != 0) {
						size++;
					}
				}
				if (c == '\n') {
					return size;
				}
				if (size == target.length) {
					return size + 1; // Indicates that target buffer is full and more stuff in the line.
				}
			}
			bufferFill = stream.read(buffer);
			pos = 0;
		}
		return size;
	}

	@Override
	public void close() throws IOException {
		if (stream != null) {
			stream.close();
		}
	}
}
