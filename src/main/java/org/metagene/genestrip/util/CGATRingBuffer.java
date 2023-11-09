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
		int a = 1;
		int b = 1;
		int count = 0;
		for (; b < maxDust; count++) {
			int h = a + b;
			a = b;
			b = h;
		}
		fib = new int[count];
		if (count > 0) {
			fib[0] = 1;
		}
		if (count > 1) {
			fib[1] = 1;
		}
		for (int i = 2; i < fib.length; i++) {
			fib[i] = fib[i - 1] + fib[i - 2];
		}
		data = new byte[size];
		seqMarks = new int[size];
		reset();
	}

	public void put(byte c) {
		data[end] = c;
		if (!(c == 'C' || c == 'G' || c == 'A' || c == 'T')) { // Inlined CGAT.isCGAT(c) for efficiency reasons.
			invalidCPos = data.length;

			sumDust = 0;
			Arrays.fill(seqMarks, 0);
			seqCount = 0;
		} else {
			if (invalidCPos > 0) {
				invalidCPos--;
			}
			if (maxDust >= 0) {
				if (c == lastChar) {
					sumDust += fib[seqCount];
					seqCount++;
					seqMarks[(end - seqCount + data.length) % data.length]++;
				} else {
					sumDust = 0;
					seqCount = 0;
				}
				lastChar = c;
			}
		}
		end = (end + 1) % data.length;
		if (maxDust >= 0) {
			int oldCount = seqMarks[(end + 1) % data.length];
			if (oldCount > 0) {
				sumDust -= fib[oldCount];
				seqMarks[(end + 2) % data.length] = oldCount - 1;
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
		sumDust = 0;
		Arrays.fill(seqMarks, 0);
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

	public void toPrintStream(PrintStream stream) {
		for (int i = 0; i < data.length; i++) {
			stream.print((char) get(i));
		}
	}
}