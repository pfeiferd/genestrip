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

public abstract class AbstractCGATBloomFilter implements Serializable {
	private static final long serialVersionUID = 1L;

	protected final static double LOG_2 = Math.log(2);

	protected static final int[] CGAT_JUMP_TABLE;
	protected static final int[] CGAT_REVERSE_JUMP_TABLE;
	
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

	protected final int k;
	protected long[] bits;
	protected long expectedInsertions;
	protected double fpp;
	protected int hashes;

	public AbstractCGATBloomFilter(int k, long expectedInsertions, double fpp) {
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
		
		hashes = optimalNumOfHashFunctions(expectedInsertions, bits);

		int logbits = (int) (Math.log((bits + 63) / 64) / LOG_2);
		int size = (1 << (logbits + 1)) - 1; // Now: size = 2^x - 1 such that 2^x > (bits + 63) / 64 < 2^(x-1)

		this.bits = new long[size];
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
	
	public int getByteSize() {
		return bits.length;
	}
	
	protected int optimalNumOfHashFunctions(long n, long m) {
		return Math.max(1, (int) Math.round((double) m / n * Math.log(2)));
	}

	protected long optimalNumOfBits(long n, double p) {
		return (long) (-n * Math.log(p) / (Math.log(2) * Math.log(2)));
	}

	public abstract void put(CGATRingBuffer buffer, int[] badPos);

	public abstract void put(byte[] seq, int start, int[] badPos);

	public abstract boolean contains(byte[] seq, int start, boolean reverse, int[] badPos);

	public abstract boolean contains(CGATRingBuffer buffer, boolean reverse, int[] badPos);
}