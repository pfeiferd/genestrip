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

import org.metagene.genestrip.util.CGATRingBuffer;

public abstract class AbstractCGATBloomFilter implements Serializable {
	public static long MAX_CAPACITY = Integer.MAX_VALUE - 8;	
	
	private static final long serialVersionUID = 1L;

	protected final static double LOG_2 = Math.log(2);

	protected final int k;
	protected final double fpp;
	protected long size;
	protected long[] bits;
	protected int hashes;
	protected long expectedInsertions;

	public AbstractCGATBloomFilter(int k, double fpp) {
		if (k <= 0) {
			throw new IllegalArgumentException("k-mer length k must be > 0");
		}
		if (fpp <= 0 || fpp >= 1) {
			throw new IllegalArgumentException("fpp must be a probability");
		}
		this.fpp = fpp;
		this.k = k;
	}

	public void clearAndEnsureCapacity(long expectedInsertions) {
		if (expectedInsertions <= 0) {
			throw new IllegalArgumentException("expected insertions must be > 0");
		}
		if (expectedInsertions > getMaxCapacity()) {
			throw new IllegalArgumentException("Expected insertions above max capacity");
		}
		this.expectedInsertions = expectedInsertions;
		
		if (size >= expectedInsertions) {
			clearArray();
		} else {
			long bits = optimalNumOfBits(expectedInsertions, fpp);
			size = (bits + 63) / 64;
			hashes = optimalNumOfHashFunctions(expectedInsertions, bits);

			initBitArray();
		}
	}
	
	public long getMaxCapacity() {
		return MAX_CAPACITY;
	}
	
	protected void clearArray() {
		Arrays.fill(bits, 0);		
	}

	protected void initBitArray() {
		this.bits = new long[(int) size];
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

	public abstract boolean put(CGATRingBuffer buffer);

	public abstract boolean put(byte[] seq, int start);

	public abstract boolean containsReverse(byte[] seq, int start, int[] badPos);
	
	public abstract boolean containsStraight(byte[] seq, int start, int[] badPos);
	
	public abstract boolean containsReverse(CGATRingBuffer buffer);

	public abstract boolean containsStraight(CGATRingBuffer buffer);		
}