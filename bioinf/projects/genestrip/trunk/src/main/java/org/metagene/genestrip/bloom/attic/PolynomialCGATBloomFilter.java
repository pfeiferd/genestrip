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
package org.metagene.genestrip.bloom.attic;

import java.util.Random;

import org.metagene.genestrip.bloom.AbstractCGATBloomFilter;
import org.metagene.genestrip.util.CGATRingBuffer;

public class PolynomialCGATBloomFilter extends AbstractCGATBloomFilter {
	private static final long serialVersionUID = 1L;
	
	private final int sizeInBits;
	private final int[][] coefficients;
	private final long p;
	
	public PolynomialCGATBloomFilter(int k, long expectedInsertions, double fpp) {
		super(k, expectedInsertions, fpp);
		sizeInBits = bits.length * 64;
		
		coefficients = new int[hashes][hashes];
		Random r = new Random(42);
		for (int i = 0; i < coefficients.length; i++) {
			int[] c = coefficients[i];
			for (int j = 0; j < coefficients.length; j++) {
				c[j] = Math.abs(r.nextInt());
			}
		}
		
		p = (1 << 31) - 1;
	}

	public void put(CGATRingBuffer buffer) {
		putViaHash(hash(buffer, false));
	}

	public void put(byte[] seq, int start) {
		putViaHash(hash(seq, start, false));
	}

	public int hash(byte[] seq, int start, boolean reverse) {
		int hash = 0;
		if (reverse) {
			for (int i = 0; i < k; i++) {
				hash = hash << 2;
				hash += CGAT_REVERSE_JUMP_TABLE[seq[k - i - 1]];
			}
		} else {
			for (int i = 0; i < k; i++) {
				hash = hash << 2;
				hash += CGAT_JUMP_TABLE[seq[i]];
			}
		}
		return hash;
	}

	public int hash(CGATRingBuffer buffer, boolean reverse) {
		assert (buffer.getSize() == k);

		int hash = 0;
		if (reverse) {
			for (int i = 0; i < k; i++) {
				hash = hash << 2;
				hash += CGAT_REVERSE_JUMP_TABLE[buffer.get(k - i - 1)];
			}
		} else {
			for (int i = 0; i < k; i++) {
				hash = hash << 2;
				hash += CGAT_JUMP_TABLE[buffer.get(i)];
			}
		}
		return hash;
	}

	private void putViaHash(int hash) {
		int index, res;
		int[] c;
		long h;
		long sum;
		
		for (int i = 0; i < coefficients.length; i++) {
			c = coefficients[i];
			sum = 0;
			h = 1;
			for (int j = 0; j < c.length; j++) {
				sum += h * c[j];
				h *= hash;
			}
			res = (int) ((sum % p) % sizeInBits);
			if (res < 0) {
				res = -res;
			}
						
			index = res >>> 6;
			bits[index] = bits[index] | (1L << (res & 0b111111));
		}
	}
	
	public boolean contains(byte[] seq, int start, boolean reverse) {
		return containsHash(hash(seq, start, reverse));
	}

	private boolean containsHash(int hash) {
		int index, res;
		int[] c;
		long h;
		long sum;
		
		for (int i = 0; i < coefficients.length; i++) {
			c = coefficients[i];
			sum = 0;
			h = 1;
			for (int j = 0; j < c.length; j++) {
				sum += h * c[j];
				h *= hash;
			}
			res = (int) ((sum % p) % sizeInBits);
			if (res < 0) {
				res = -res;
			}
						
			index = res >>> 6;
			if (((bits[index] >> (res & 0b111111)) & 1L) == 0) {
				return false;
			}
		}
		return true;
	}

	public boolean contains(CGATRingBuffer buffer, boolean reverse) {
		return containsHash(hash(buffer, reverse));
	}
}