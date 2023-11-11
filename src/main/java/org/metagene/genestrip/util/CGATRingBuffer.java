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
package org.metagene.genestrip.util;

import java.io.PrintStream;
import java.io.Serializable;
import java.util.Arrays;

public class CGATRingBuffer implements Serializable {
	private static final long serialVersionUID = 1L;

	private final int seqMarks[];
	private final byte[] data;

	private int end;
	private boolean filled;
	private int invalidCPos;

	private final int maxDust;
	private final int[] fib;

	private byte lastChar;
	private int sumDust;
	private int seqCount;

	public CGATRingBuffer(int size) {
		this(size, -1);
	}

	public CGATRingBuffer(int size, int maxDust) {
		if (maxDust > Short.MAX_VALUE) {
			throw new IllegalArgumentException("Unreasonable size for maxDust :" + maxDust);
		}
		this.maxDust = maxDust;
		if (maxDust >= 0) {
			fib = new int[size];
			if (size > 0) {
				fib[0] = 1;
			}
			if (size > 1) {
				fib[1] = 1;
			}
			if (size > 2) {
				fib[2] = 1;
			}
			for (int i = 3; i < fib.length; i++) {
				fib[i] = fib[i - 1] + fib[i - 2];
			}
			seqMarks = new int[size];
		}
		else {
			seqMarks = null;
			fib = null;
		}
		data = new byte[size];
		reset();
	}

	public void put(byte c) {
		data[end] = c;
		if (!(c == 'C' || c == 'G' || c == 'A' || c == 'T')) { // Inlined CGAT.isCGAT(c) for efficiency reasons.
			invalidCPos = data.length;
			sumDust = maxDust >= 0 ? 0 : -1;
		} else {
			if (invalidCPos > 0) {
				invalidCPos--;
			}
			if (maxDust >= 0) {
				if (c == lastChar) {
					seqMarks[(end - seqCount + data.length) % data.length]++;
					seqCount++;
					sumDust += fib[seqCount - 1];
				} else {
					seqCount = 0;
				}
				lastChar = c;
			}
		}
		end = (end + 1) % data.length;
		if (maxDust >= 0) {
			int oldCount = seqMarks[end];
			seqMarks[end % data.length] = 0;
			if (oldCount > 0) {
				sumDust -= fib[oldCount - 1];
				seqMarks[(end + 1) % data.length] = oldCount - 1;
			}
		}
		if (end == 0) {
			filled = true;
		}
	}

	public void reset() {
		end = 0;
		invalidCPos = 0;
		filled = false;

		seqCount = 0;
		sumDust = maxDust >= 0 ? 0 : -1;
		if (seqMarks != null) {
			Arrays.fill(seqMarks, 0);
		}
	}

	public boolean isFilled() {
		return filled;
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

	public boolean isDust() {
		return sumDust > maxDust;
	}

	public int getDustValue() {
		return sumDust;
	}
	
	public int getMaxDust() {
		return maxDust;
	}

	public void toPrintStream(PrintStream stream) {
		for (int i = 0; i < data.length; i++) {
			stream.print((char) get(i));
		}
	}
}