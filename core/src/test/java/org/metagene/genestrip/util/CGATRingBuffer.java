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

public class CGATRingBuffer extends CGATLongBuffer {
	private final byte[] data;

	public CGATRingBuffer(int size, int maxDust) {
		super(size, maxDust);
		data = new byte[size];
	}

	public CGATRingBuffer(int size) {
		this(size, -1);
	}

	public long putForTest(byte c) {
		data[bpCounter] = c;
		return super.put(c);
	}

	public byte get(int index) {
		return data[(bpCounter + index) % data.length];
	}
	
	public long getKMerForTest() {
		if (!filled) {
			return -1;
		}
		
		long res = 0;
		int c;
		for (int i = 0; i < size; i++) {
			// Inlined: res = Long.rotateLeft(res, 2);
			res = (res << 2) | (res >>> -2);
			c = CGAT.CGAT_JUMP_TABLE[data[(bpCounter + i) % size]];
			res += c;
		}

		return res;
	}
	
	public String toString() {
		StringBuilder builder = new StringBuilder();

		for (int i = 0; i < size; i++) {
			builder.append((char)data[(bpCounter + i) % size]);
		}
		/* For debug purposes
		if (srl0Buffer != null) {
			builder.append(" [");
			for (int i = 0; i < size; i++) {
				if (i > 0) {
					builder.append(',');
				}
				builder.append(srl1Buffer[(bpCounter + i) % size]);
			}
			builder.append(']');
		}
		 */

		return builder.toString();
	}
}
