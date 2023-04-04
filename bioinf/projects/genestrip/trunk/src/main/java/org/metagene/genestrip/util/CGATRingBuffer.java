package org.metagene.genestrip.util;

import java.io.PrintStream;
import java.io.Serializable;

public class CGATRingBuffer implements Serializable {
	private static final long serialVersionUID = 1L;

	// Made public for fast access:
	public int end;
	public byte[] data;
	public byte[] directPut;
	public int directPutStart;
	public boolean filled;

	private int invalidCPos;

	public CGATRingBuffer(int size) {
		data = new byte[size];
		reset();
	}

	public void put(byte c) {
		data[end] = c;
		if (isInvalidByte(c)) {
			invalidCPos = data.length;
		} else if (invalidCPos > 0) {
			invalidCPos--;
		}
		end = (end + 1) % data.length;
		if (end == 0) {
			filled = true;
		}
	}
	
	public void reset() {
		end = 0;
		invalidCPos = 0;
		filled = false;
	}
	
	public boolean isFilled() {
		return filled;
	}
	
	protected boolean isInvalidByte(byte b) {
		return CGAT.isCGAT(b);
	}

	public int getSize() {
		return data.length;
	}

	public byte get(int index) {
		return data[(end + index) % data.length];
	}

	public String toString() {
		StringBuilder builder = new StringBuilder();

		for (int i = 0; i < data.length; i++) {
			builder.append((char) get(i));
		}

		return builder.toString();
	}

	public boolean isCGAT() {
		return invalidCPos == 0;
	}

	public void toStream(PrintStream stream) {
		for (int i = 0; i < data.length; i++) {
			stream.append((char) get(i));
		}
	}
}