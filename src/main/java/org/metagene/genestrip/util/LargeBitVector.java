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

import java.util.Arrays;

import it.unimi.dsi.fastutil.BigArrays;

public class LargeBitVector {
	public static long MAX_SMALL_CAPACITY = Integer.MAX_VALUE - 8;

	private final long size;
	protected final long[] bits;
	protected final long[][] largeBits;

	public LargeBitVector(long size) {
		this.size = (size + 63) / 64;
		if (size > MAX_SMALL_CAPACITY) {
			bits = null;
			largeBits = BigArrays.ensureCapacity(BigArrays.wrap(new long[0]), size);
		} else {
			bits = new long[(int) size];
			largeBits = null;
		}
	}
	
	public void clear() {
		if (largeBits != null) {
			BigArrays.fill(largeBits, 0);
		} else {
			Arrays.fill(bits, 0);
		}		
	}
	
	public long getSizeBitSize() {
		return size * 64;
	}

	public void set(long index) {
		if (largeBits != null) {
			long arrayIndex = ((index >>> 6) % size);
			BigArrays.set(largeBits, arrayIndex, BigArrays.get(largeBits, arrayIndex) | (1L << (index & 0b111111)));
		} else {
			int arrayIndex = (int) ((index >>> 6) % size);
			bits[arrayIndex] = bits[arrayIndex] | (1L << (index & 0b111111));
		}
	}

	public boolean get(long index) {
		if (largeBits != null) {
			long arrayIndex = ((index >>> 6) % size);
			return ((BigArrays.get(largeBits, arrayIndex) >> (index & 0b111111)) & 1L) == 1;
		} else {
			int arrayIndex = (int) ((index >>> 6) % size);
			return (((bits[arrayIndex] >> (index & 0b111111)) & 1L) == 1);
		}
	}
}
