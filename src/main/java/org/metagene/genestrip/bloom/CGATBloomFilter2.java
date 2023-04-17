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
import java.util.Random;

import org.metagene.genestrip.util.CGATRingBuffer;

public class CGATBloomFilter2 implements Serializable {
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
	private final long[] bits;
	private final int sizeInBits;
	private final long expectedInsertions;
	private final double fpp;
	
	// For tabulation hashing
	private final int[][][] t;
	private final int permBytes;
	private final long permMax;
	private final long p;
	
	public CGATBloomFilter2(int k, long expectedInsertions, double fpp) {
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
		int size = (int) ((bits + 63L) / 64L); // Size as number of longs (which has 64 bits);
		this.bits = new long[size];
		sizeInBits = size * 64;
		int hashFunctions = optimalNumOfHashFunctions(expectedInsertions, sizeInBits);
		
		permBytes = (int) ((Math.log(size) / LOG_2) + 7) / 8;
		t = new int[hashFunctions][permBytes][256];
		Random r = new Random(42);
		for (int i = 0; i < t.length; i++) {
			int[][] perm = t[i];
			for (int j = 0; j < perm.length; j++) {
				for (int l = 0; l < 256; l++) {
					perm[j][l] = r.nextInt();					
				}
			}
		}
		permMax = 1L << (8 * permBytes);
		
		p = ((permMax << 3) - 1); // Larger and relatve prime should suffice.
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
		return t.length;
	}

	private int optimalNumOfHashFunctions(long n, long m) {
		return Math.max(1, (int) Math.round((m * LOG_2) / n));
	}

	private long optimalNumOfBits(long n, double p) {
		return (long) (-n * Math.log(p) / (LOG_2 * LOG_2));
	}

	public void put(CGATRingBuffer buffer) {
		long hash = hash(buffer, false);
		putViaHash(hash);
	}

	public void put(byte[] seq, int start) {
		long hash = hash(seq, start, false);
		putViaHash(hash);
	}

	public long hash(byte[] seq, int start, boolean reverse) {
		long hash = 0;
		if (reverse) {
			for (int i = 0; i < k; i++) {
				hash += CGAT_REVERSE_JUMP_TABLE[seq[k - i - 1]];
				hash = hash * 3;
			}
		} else {
			for (int i = 0; i < k; i++) {
				hash += CGAT_JUMP_TABLE[seq[i]];
				hash = hash * 3;
			}
		}
		return hash;
	}

	public long hash(CGATRingBuffer buffer, boolean reverse) {
		assert (buffer.getSize() == k);

		long hash = 0;
		if (reverse) {
			for (int i = 0; i < k; i++) {
				hash += CGAT_REVERSE_JUMP_TABLE[buffer.get(k - i - 1)];
				hash = hash * 3;
			}
		} else {
			for (int i = 0; i < k; i++) {
				hash += CGAT_JUMP_TABLE[buffer.get(i)];
				hash = hash * 3;
			}
		}
		return hash;
	}

	private void putViaHash(long hash) {
		int nextHash;
		int index;
		int[][] perm;
		
		int reducedHash = (int) ((hash % p) % permMax);
		if (reducedHash < 0) {
			reducedHash = -reducedHash;
		}
		int r;
		
		for (int i = 0; i < t.length; i++) {
			perm = t[i];
			nextHash = 0;
			r = reducedHash;
			for (int j = 0; j < perm.length; j++) {
				nextHash += perm[j][r & 0b11111111];
				r >>= 8;
			}
			nextHash = nextHash % sizeInBits;
			if (nextHash < 0) {
				nextHash = -nextHash;
			}
						
			index = nextHash >>> 6;
			bits[index] = bits[index] | (1L << (nextHash & 0b111111));
		}
	}
	
	public boolean contains(byte[] seq, int start, boolean reverse) {
		long hash = hash(seq, start, reverse);
		return containsHash(hash);
	}

	private boolean containsHash(long hash) {
		int nextHash;
		int index;
		int[][] perm;
		
		int reducedHash = (int) ((hash % p) % permMax);
		if (reducedHash < 0) {
			reducedHash = -reducedHash;
		}
		int r;
		
		for (int i = 0; i < t.length; i++) {
			perm = t[i];
			nextHash = 0;
			r = reducedHash;
			for (int j = 0; j < perm.length; j++) {
				nextHash += perm[j][r & 0b11111111];
				r >>= 8;
			}
			nextHash = nextHash % sizeInBits;
			if (nextHash < 0) {
				nextHash = -nextHash;
			}
						
			index = nextHash >>> 6;
			if (((bits[index] >>> (nextHash & 0b111111)) & 1L) == 0) {
				return false;
			}
		}
		return true;
	}

	public boolean contains(CGATRingBuffer buffer, boolean reverse) {
		long hash = hash(buffer, reverse);
		return containsHash(hash);
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