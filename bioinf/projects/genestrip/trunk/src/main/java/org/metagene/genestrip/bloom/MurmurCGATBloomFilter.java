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
package org.metagene.genestrip.bloom;

import java.io.Serializable;
import java.util.Arrays;
import java.util.Random;

import org.metagene.genestrip.util.CGAT;
import org.metagene.genestrip.util.CGATRingBuffer;

import it.unimi.dsi.fastutil.BigArrays;

public class MurmurCGATBloomFilter implements Serializable {
	public static long MAX_SMALL_CAPACITY = Integer.MAX_VALUE - 8;

	private static final long serialVersionUID = 1L;

	protected final int k;
	protected final double fpp;
	protected final Random random;

	protected long expectedInsertions;
	protected long size;

	protected long[] bits;
	protected long[][] largeBits;
	protected boolean large;

	protected int hashes;
	protected long[] hashFactors;

	public MurmurCGATBloomFilter(int k, double fpp) {
		if (k <= 0) {
			throw new IllegalArgumentException("k-mer length k must be > 0");
		}
		if (fpp <= 0 || fpp >= 1) {
			throw new IllegalArgumentException("fpp must be a probability");
		}
		this.fpp = fpp;
		this.k = k;

		large = false;
		random = new Random(42);
	}

	public void clearAndEnsureCapacity(long expectedInsertions) {
		if (expectedInsertions <= 0) {
			throw new IllegalArgumentException("expected insertions must be > 0");
		}
		this.expectedInsertions = expectedInsertions;

		long bits = optimalNumOfBits(expectedInsertions, fpp);
		long newSize = (bits + 63) / 64;
		if (size >= newSize) {
			clearArray();
		} else {
			size = newSize;
			hashes = optimalNumOfHashFunctions(expectedInsertions, bits);
			hashFactors = new long[hashes];
			for (int i = 0; i < hashFactors.length; i++) {
				hashFactors[i] = random.nextLong();
			}

			initBitArray();
		}
	}

	protected void initBitArray() {
		if (size > MAX_SMALL_CAPACITY || large == true) {
			large = true;
			bits = null;
			if (largeBits == null) {
				largeBits = BigArrays.ensureCapacity(BigArrays.wrap(new long[0]), size);
			} else {
				largeBits = BigArrays.ensureCapacity(largeBits, size);
			}
		} else {
			this.bits = new long[(int) size];
		}
	}

	protected void clearArray() {
		if (large) {
			BigArrays.fill(largeBits, 0);
		} else {
			Arrays.fill(bits, 0);
		}
	}

	public long getExpectedInsertions() {
		return expectedInsertions;
	}

	public int getK() {
		return k;
	}

	public int getHashes() {
		return hashes;
	}

	public double getFpp() {
		return fpp;
	}

	public long getByteSize() {
		return size * 8;
	}

	protected int optimalNumOfHashFunctions(long n, long m) {
		return Math.max(1, (int) Math.round((double) m / n * Math.log(2)));
	}

	protected long optimalNumOfBits(long n, double p) {
		return (long) (-n * Math.log(p) / (Math.log(2) * Math.log(2)));
	}

	public boolean isLarge() {
		return large;
	}

	public void makeLarge() {
		this.large = true;
	}

	public void put(CGATRingBuffer buffer) {
		putViaHash(CGAT.kmerToLongStraight(buffer));
	}

	public void put(byte[] seq, int start) {
		putViaHash(CGAT.kMerToLongStraight(seq, start, k, null));
	}

	public void putLong(long data) {
		putViaHash(data);
	}
	
	protected void putViaHash(long data) {
		int index;
		long hash;

		if (large) {
			for (int i = 0; i < hashes; i++) {
				hash = hash(data, i);
				index = (int) ((hash >>> 6) % bits.length);
				BigArrays.set(largeBits, index, BigArrays.get(largeBits, index) | (1L << (hash & 0b111111)));
			}
		} else {
			for (int i = 0; i < hashes; i++) {
				hash = hash(data, i);
				index = (int) ((hash >>> 6) % bits.length);
				bits[index] = bits[index] | (1L << (hash & 0b111111));
			}
		}
	}

	public boolean containsStraight(byte[] seq, int start, int[] badPos) {
		long data = CGAT.kMerToLongStraight(seq, start, k, badPos);
		if (data == -1 && badPos != null && badPos[0] == -1) {
			return false;
		}
		return containsViaHash(data);
	}

	public boolean containsReverse(byte[] seq, int start, int[] badPos) {
		long data = CGAT.kMerToLongReverse(seq, start, k, badPos);
		if (data == -1 && badPos != null && badPos[0] == -1) {
			return false;
		}
		return containsViaHash(data);
	}

	protected boolean containsViaHash(long data) {
		long hash;
		int index;

		if (large) {
			for (int i = 0; i < hashes; i++) {
				hash = hash(data, i);
				index = (int) ((hash >>> 6) % bits.length);
				if (((BigArrays.get(largeBits, index) >> (hash & 0b111111)) & 1L) == 0) {
					return false;
				}
			}
		} else {
			for (int i = 0; i < hashes; i++) {
				hash = hash(data, i);
				index = (int) ((hash >>> 6) % bits.length);
				if (((bits[index] >> (hash & 0b111111)) & 1L) == 0) {
					return false;
				}
			}
		}
		return true;
	}

	public boolean containsStraight(CGATRingBuffer buffer) {
		long data = CGAT.kmerToLongStraight(buffer);
		if (data == 0) {
			return false;
		}
		return containsViaHash(data);
	}

	public boolean containsReverse(CGATRingBuffer buffer) {
		long data = CGAT.kmerToLongReverse(buffer);
		if (data == 0) {
			return false;
		}
		return containsViaHash(data);
	}
	
	public boolean containsLong(long data) {
		if (data == 0) {
			return false;
		}
		return containsViaHash(data);		
	}

	/*
	 * The following code was copied from
	 * org.apache.commons.codec.digest.MurmurHash3.hash64()
	 * 
	 * The original function is deprecated and does not support changing seed.
	 * That's why its implemented here...
	 */

	// Constants for 128-bit variant
	private static final long C1 = 0x87c37b91114253d5L;
	private static final long C2 = 0x4cf5ad432745937fL;
	private static final int R1 = 31;
	private static final int R2 = 27;
	private static final int M = 5;
	private static final int N1 = 0x52dce729;

	protected long hash(long data, int i) {
		long hash = hashFactors[i];

		// Reverse byte inlined from Long.reverseBytes()
		long k = (data & 0x00ff00ff00ff00ffL) << 8 | (data >>> 8) & 0x00ff00ff00ff00ffL;
		k = (k << 48) | ((k & 0xffff0000L) << 16) | ((k >>> 16) & 0xffff0000L) | (k >>> 48);

		final int length = Long.BYTES;
		// mix functions
		k *= C1;
		// Rotate left inlined from Long.rotateLeft()
		k = (k << R1) | (k >>> -R1);
		k *= C2;
		hash ^= k;
		// Rotate left inlined from Long.rotateLeft()
		hash = ((hash << R2) | (hash >>> -R2)) * M + N1;
		// finalization
		hash ^= length;
		
		// Inlined from MurmurHash3.fmix64()
		hash ^= (hash >>> 33);
		hash *= 0xff51afd7ed558ccdL;
		hash ^= (hash >>> 33);
		hash *= 0xc4ceb9fe1a85ec53L;
		hash ^= (hash >>> 33);
		return hash;
	}
}