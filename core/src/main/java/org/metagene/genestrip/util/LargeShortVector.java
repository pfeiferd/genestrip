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

/**
 * A vector of short counters that transparently switches between a single {@code short[]} for small
 * capacities and fastutil big arrays for capacities beyond the range of a normal array, and never
 * shrinks.
 */
public class LargeShortVector implements Serializable {
	private static final long serialVersionUID = 1L;

	/** Maximum capacity that still uses the small ({@code int}-indexed) array. */
	public static long MAX_SMALL_CAPACITY = Integer.MAX_VALUE - 8;

	/** Current size (number of shorts) of the vector. */
	protected long size;
	// Made public for inlining
	/** Small ({@code int}-indexed) backing array, or {@code null} when large storage is used. */
	public short[] shorts;
	// Made public for inlining
	/** Large ({@code long}-indexed) backing array, or {@code null} when small storage is used. */
	public short[][] largeShorts;

	/**
	 * Creates a vector with at least the given initial capacity (number of shorts).
	 *
	 * @param initialSize the initial capacity in shorts
	 */
	public LargeShortVector(long initialSize) {
		this(initialSize, false);
	}

	/**
	 * Creates a vector with at least the given initial capacity (number of shorts).
	 *
	 * @param initialSize  the initial capacity in shorts
	 * @param enforceLarge if true, forces the big-array backing even when a normal array would fit.
	 */
	public LargeShortVector(long initialSize, boolean enforceLarge) {
		size = 0;
		ensureCapacity(initialSize, enforceLarge);
	}

	/**
	 * Sets all shorts to zero.
	 */
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
	 * @param enforceLarge if true, forces the big-array backing even when a normal array would fit.
	 * @return Whether the vector got bigger or not.
	 */
	public boolean ensureCapacity(long newSize, boolean enforceLarge) {
		if (newSize > size) {
			size = newSize;
			if (size > MAX_SMALL_CAPACITY || enforceLarge || largeShorts != null) {
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

	/**
	 * Returns whether this vector currently uses the large (big-array) backing storage.
	 *
	 * @return {@code true} if the large backing storage is in use
	 */
	public boolean isLarge() {
		return largeShorts != null;
	}

	/**
	 * Returns the current size (number of shorts) of the vector.
	 *
	 * @return the current size in shorts
	 */
	public long getSize() {
		return size;
	}
	
	/**
	 * Increments the short at the given index and returns its new value.
	 *
	 * @param index the index of the short to increment
	 * @return the incremented value
	 */
	public final short inc(final long index) {
		if (largeShorts != null) {
			BigArrays.incr(largeShorts, index);
			return BigArrays.get(largeShorts, index);
		} else {
			return ++shorts[(int) index];
		}
	}

	/**
	 * Returns the short at the given index.
	 *
	 * @param index the index of the short to read
	 * @return the short value at the index
	 */
	public final short get(long index) {
		if (largeShorts != null) {
			return BigArrays.get(largeShorts, index);
		} else {
			return shorts[(int) index];
		}
	}
}
