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
		if (!(c == 'C' || c == 'G' || c == 'A' || c == 'T')) { // Inlined CGAT.isCGAT(c) for efficiency reasons.
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

	public void toPrintStream(PrintStream stream) {
		for (int i = 0; i < data.length; i++) {
			stream.print((char) get(i));
		}
	}
}