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

/**
 * A vector of short counters that transparently switches between a single {@code short[]} for small
 * capacities and a self-managed array of {@code short[]} buckets for capacities beyond the range of a
 * normal array, and never shrinks.
 * <p>
 * The large backing is a plain {@code short[][]} (an "array of arrays") laid out on a fixed grid:
 * every bucket holds exactly {@link #BUCKET_SIZE} shorts except the last, which holds the remainder.
 * The grid is self-managed (no external big-array library is involved); the bucket width matches the
 * historic segmentation, so the serialized form is unchanged.
 */
public class LargeShortVector implements Serializable {
	private static final long serialVersionUID = 1L;

	/** Maximum capacity that still uses the small ({@code int}-indexed) array. */
	public static long MAX_SMALL_CAPACITY = Integer.MAX_VALUE - 8;

	/** Base-2 logarithm of {@link #BUCKET_SIZE}; an index is split into bucket and displacement here. */
	private static final int BUCKET_SHIFT = 27;
	/** Number of shorts per bucket of the large backing (all but the last are full). */
	private static final int BUCKET_SIZE = 1 << BUCKET_SHIFT;
	/** Mask selecting the in-bucket displacement of an index. */
	private static final int BUCKET_MASK = BUCKET_SIZE - 1;

	/** Current size (number of shorts) of the vector. */
	protected long size;
	// Made public for inlining
	/** Small ({@code int}-indexed) backing array, or {@code null} when large storage is used. */
	public short[] shorts;
	// Made public for inlining
	/** Large (bucketed) backing array, or {@code null} when small storage is used. */
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
			for (short[] bucket : largeShorts) {
				Arrays.fill(bucket, (short) 0);
			}
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
			long oldCount = size;
			size = newSize;
			if (size > MAX_SMALL_CAPACITY || enforceLarge || largeShorts != null) {
				largeShorts = growLarge(shorts, largeShorts, oldCount, size);
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
	 * Allocates the large (bucketed) backing for {@code newCount} shorts on the fixed {@link
	 * #BUCKET_SIZE} grid and copies the first {@code oldCount} shorts of the previous backing into it.
	 * Exactly one of {@code smallOld} / {@code largeOld} is non-null (the previous backing), or both
	 * are null on the very first allocation.
	 *
	 * @param smallOld the previous small backing, or {@code null}
	 * @param largeOld the previous large backing, or {@code null}
	 * @param oldCount the number of shorts to preserve from the previous backing
	 * @param newCount the desired capacity in shorts
	 * @return the freshly allocated large backing holding the preserved shorts
	 */
	private static short[][] growLarge(short[] smallOld, short[][] largeOld, long oldCount, long newCount) {
		int bucketCount = (int) ((newCount + BUCKET_SIZE - 1) >>> BUCKET_SHIFT);
		short[][] result = new short[bucketCount][];
		for (int b = 0; b < bucketCount; b++) {
			long start = (long) b << BUCKET_SHIFT;
			result[b] = new short[(int) Math.min(BUCKET_SIZE, newCount - start)];
		}
		// Copy the retained shorts bucket-aligned. Source and destination share the same grid, so a
		// block never straddles more than one source and one destination bucket.
		long copied = 0;
		while (copied < oldCount) {
			int dstBucket = (int) (copied >>> BUCKET_SHIFT);
			int dstOff = (int) (copied & BUCKET_MASK);
			short[] src;
			int srcOff;
			int srcRoom;
			if (largeOld != null) {
				src = largeOld[(int) (copied >>> BUCKET_SHIFT)];
				srcOff = (int) (copied & BUCKET_MASK);
				srcRoom = src.length - srcOff;
			} else {
				src = smallOld;
				srcOff = (int) copied;
				srcRoom = smallOld.length - srcOff;
			}
			int n = (int) Math.min(Math.min(result[dstBucket].length - dstOff, srcRoom), oldCount - copied);
			System.arraycopy(src, srcOff, result[dstBucket], dstOff, n);
			copied += n;
		}
		return result;
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
			return ++largeShorts[(int) (index >>> BUCKET_SHIFT)][(int) (index & BUCKET_MASK)];
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
			return largeShorts[(int) (index >>> BUCKET_SHIFT)][(int) (index & BUCKET_MASK)];
		} else {
			return shorts[(int) index];
		}
	}
}
