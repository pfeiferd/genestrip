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

import java.io.Serializable;
import java.util.Arrays;

import it.unimi.dsi.fastutil.BigArrays;
import it.unimi.dsi.fastutil.shorts.ShortBigArrays;

public class LargeShortVector implements Serializable {
	private static final long serialVersionUID = 1L;

	public static long MAX_SMALL_CAPACITY = Integer.MAX_VALUE - 8;

	protected long size;
	// Made public for inlining
	public short[] shorts;
	// Made public for inlining
	public short[][] largeShorts;

	public LargeShortVector(long initialSize) {
		this(initialSize, false);
	}

	public LargeShortVector(long initialSize, boolean enforceLarge) {
		size = 0;
		ensureCapacity(initialSize, enforceLarge);
	}

	public void clear() {
		if (largeShorts != null) {
			BigArrays.fill(largeShorts, (short) 0);
		} else if (shorts != null) {
			Arrays.fill(shorts, (short) 0);
		}
	}

	/**
	 * Ensure that the bit vector has at least the desired size. The bit vector
	 * might be enlarged but never gets smaller in size.
	 * 
	 * @param newSize The desired size in bits.
	 * @return Whether the vector got bigger or not.
	 */
	public boolean ensureCapacity(long newSize, boolean enforceLarge) {
		if (newSize > size) {
			size = newSize;
			if (size > MAX_SMALL_CAPACITY || enforceLarge) {
				if (shorts != null) {
					largeShorts = BigArrays.wrap(shorts);
				}
				largeShorts = BigArrays
						.ensureCapacity(largeShorts == null ? ShortBigArrays.EMPTY_BIG_ARRAY : largeShorts, size);
				shorts = null;
			} else {
				short[] oldShorts = shorts;
				shorts = new short[(int) size];
				if (oldShorts != null) {
					System.arraycopy(oldShorts, 0, shorts, 0, oldShorts.length);
				}
				largeShorts = null;
			}
			return true;
		} else {
			return false;
		}
	}

	public boolean isLarge() {
		return largeShorts != null;
	}

	public long getSize() {
		return size;
	}
	
	public final short inc(final long index) {
		if (largeShorts != null) {
			BigArrays.incr(largeShorts, index);
			return BigArrays.get(largeShorts, index);
		} else {
			return ++shorts[(int) index];
		}
	}

	public final short get(long index) {
		if (largeShorts != null) {
			return BigArrays.get(largeShorts, index);
		} else {
			return shorts[(int) index];
		}
	}
}
