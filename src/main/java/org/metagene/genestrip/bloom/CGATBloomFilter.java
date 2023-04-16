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

import org.metagene.genestrip.util.CGATRingBuffer;

public class CGATBloomFilter implements Serializable {
	private static final long serialVersionUID = 1L;

	private final static double LOG_2 = Math.log(2);

	private static final int[] CGAT_JUMP_TABLE;
	private static final int[] CGAT_REVERSE_JUMP_TABLE;
	static {
		CGAT_JUMP_TABLE = new int[Byte.MAX_VALUE];
		CGAT_REVERSE_JUMP_TABLE = new int[Byte.MAX_VALUE];
		for (int i = 0; i < CGAT_JUMP_TABLE.length; i++) {
			CGAT_JUMP_TABLE[i] = -1;
			CGAT_REVERSE_JUMP_TABLE[i] = -1;
		}
		CGAT_JUMP_TABLE['C'] = 0;
		CGAT_JUMP_TABLE['G'] = 1;
		CGAT_JUMP_TABLE['A'] = 2;
		CGAT_JUMP_TABLE['T'] = 3;

		CGAT_REVERSE_JUMP_TABLE['C'] = 1;
		CGAT_REVERSE_JUMP_TABLE['G'] = 0;
		CGAT_REVERSE_JUMP_TABLE['A'] = 3;
		CGAT_REVERSE_JUMP_TABLE['T'] = 2;
	}

	private final int k;
	private long[] bits;
	private final int[] hashFactor;
	private long expectedInsertions;
	private double fpp;

	public CGATBloomFilter(int k, long expectedInsertions, double fpp) {
		if (k <= 0) {
			throw new IllegalArgumentException("k-mer length k must be > 0");
		}
		if (expectedInsertions <= 0) {
			throw new IllegalArgumentException("expected insertions must be > 0");
		}
		if (fpp <= 0 || fpp >= 1) {
			throw new IllegalArgumentException("fpp must be a probability");
		}
		this.fpp = fpp;
		this.k = k;

		this.expectedInsertions = expectedInsertions;
		long bits = optimalNumOfBits(expectedInsertions, fpp);

		int logbits = (int) (Math.log((bits + 63) / 64) / LOG_2);
		int size = (1 << (logbits + 1)) - 1; // Now: size = 2^x - 1 such that 2^x > (bits + 63) / 64 < 2^(x-1)

		this.bits = new long[size];
		this.hashFactor = primeNumbersBruteForce(optimalNumOfHashFunctions(expectedInsertions, size * 64));
	}
	
	public long getExpectedInsertions() {
		return expectedInsertions;
	}
	
	public int getK() {
		return k;
	}
	
	public double getFpp() {
		return fpp;
	}
	
	public int getByteSize() {
		return bits.length;
	}
	
	public int getHashFunctions() {
		return hashFactor.length;
	}

	private int optimalNumOfHashFunctions(long n, long m) {
		// (m / n) * log(2), but avoid truncation due to division!
		return Math.max(1, (int) Math.round((double) m / n * Math.log(2)));
	}

	private long optimalNumOfBits(long n, double p) {
		return (long) (-n * Math.log(p) / (Math.log(2) * Math.log(2)));
	}

	public void put(CGATRingBuffer buffer) {
		long hash1 = hash(buffer, true, false);
		long hash2 = hash(buffer, false, false);
		putViaHash(hash1, hash2);
	}

	public void put(byte[] seq, int start) {
		long hash1 = hash(seq, start, true, false);
		long hash2 = hash(seq, start, false, false);
		putViaHash(hash1, hash2);
	}

	public long hash(byte[] seq, int start, boolean even, boolean reverse) {
		long hash = 0;
		if (reverse) {
			for (int i = even ? 0 : 1; i < k; i += 2) {
				hash = hash << 2;
				hash += CGAT_REVERSE_JUMP_TABLE[seq[k - i - 1]];
			}
		} else {
			for (int i = even ? 0 : 1; i < k; i += 2) {
				hash = hash << 2;
				hash += CGAT_JUMP_TABLE[seq[i]];
			}
		}
		return hash;
	}

	public long hash(CGATRingBuffer buffer, boolean even, boolean reverse) {
		assert (buffer.getSize() == k);

		long hash = 0;
		if (reverse) {
			for (int i = even ? 0 : 1; i < k; i += 2) {
				hash = hash << 2;
				hash += CGAT_REVERSE_JUMP_TABLE[buffer.get(k - i - 1)];
			}
		} else {
			for (int i = even ? 0 : 1; i < k; i += 2) {
				hash = hash << 2;
				hash += CGAT_JUMP_TABLE[buffer.get(i)];
			}
		}
		return hash;
	}

	private void putViaHash(long hash1, long hash2) {
		long hash;
		int index;
		
		for (int i = 0; i < hashFactor.length; i++) {
			hash = hash1 * hashFactor[i] + hash2;
			index = (int) ((hash >>> 6) % bits.length);
			bits[index] = bits[index] | (1L << (hash & 0b111111L));
		}
	}

	public boolean contains(byte[] seq, int start, boolean reverse) {
		long hash1 = hash(seq, start, true, reverse);
		long hash2 = hash(seq, start, false, reverse);
		return containsHash(hash1, hash2);
	}

	private boolean containsHash(long hash1, long hash2) {
		long hash;
		int index;
		
		for (int i = 0; i < hashFactor.length; i++) {
			hash = hash1 * hashFactor[i] + hash2;
			index = (int) ((hash >>> 6) % bits.length);
			if (((bits[index] >> (hash & 0b111111L)) & 1L) == 0) {
				return false;
			}
		}
		return true;
	}

	public boolean contains(CGATRingBuffer buffer, boolean reverse) {
		long hash1 = hash(buffer, true, reverse);
		long hash2 = hash(buffer, false, reverse);
		return containsHash(hash1, hash2);
	}

	// 2 is not in it...
	public static int[] primeNumbersBruteForce(int n) {
		int[] res = new int[n];
		int count = 0;
		int num = 3;
		while(count < n) { 
			boolean prime = true;// to determine whether the number is prime or not
			int max = (int) Math.sqrt(num);
			for (int i = 2; i <= max; i++) {
				if (num % i == 0) {
					prime = false;
					break;
				}				
			}
			if (prime) {
				res[count++] = num;
			}
			num++;
		}
		return res;
	}
}