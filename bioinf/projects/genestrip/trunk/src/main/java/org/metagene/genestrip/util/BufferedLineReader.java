package org.metagene.genestrip.util;

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

	public int nextLine(byte[] target) throws IOException {
		int size = 0;
		while (bufferFill != -1) {
			if (pos < bufferFill) {
				byte c = -1;
				for (size = 0; size < target.length && pos < bufferFill && c != '\n'; size++, pos++) {
					target[size] = c = buffer[pos];
				}
				if (size == target.length || c == '\n') {
					return size;
				}
			}
			if (pos == bufferFill && bufferFill < buffer.length) {
				return size;
			}
			bufferFill = stream.read(buffer);
			pos = 0;
		}
		return size;
	}
}
