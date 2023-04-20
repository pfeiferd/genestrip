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

public class TabulationCGATBloomFilter extends AbstractCGATBloomFilter {
	private static final long serialVersionUID = 1L;
	
	// For tabulation hashing
	private final int[][][] t;
	private final int permBytes;
	private final int sizeInBits;
	
	public TabulationCGATBloomFilter(int k, long expectedInsertions, double fpp) {
		super(k, expectedInsertions, fpp);

		sizeInBits = bits.length * 64;

		permBytes = (int) ((Math.log(sizeInBits) / LOG_2) + 7) / 8;
		t = new int[hashes][permBytes][256];
		Random r = new Random(42);
		for (int i = 0; i < t.length; i++) {
			int[][] perm = t[i];
			for (int j = 0; j < perm.length; j++) {
				for (int l = 0; l < 256; l++) {
					perm[j][l] = r.nextInt();					
				}
			}
		}
	}

	public void put(CGATRingBuffer buffer, int[] badPos) {
		long hash = hash(buffer, false);
		putViaHash(hash);
	}

	public void put(byte[] seq, int start, int[] badPos) {
		long hash = hash(seq, start, false);
		putViaHash(hash);
	}

	public long hash(byte[] seq, int start, boolean reverse) {
		long hash = 0;
		if (reverse) {
			for (int i = 0; i < k; i++) {
				hash = hash << 2;
				hash ^= CGAT_REVERSE_JUMP_TABLE[seq[k - i - 1]];
			}
		} else {
			for (int i = 0; i < k; i++) {
				hash = hash << 2;
				hash ^= CGAT_JUMP_TABLE[seq[i]];
			}
		}
		return hash;
	}

	public long hash(CGATRingBuffer buffer, boolean reverse) {
		assert (buffer.getSize() == k);

		long hash = 0;
		if (reverse) {
			for (int i = 0; i < k; i++) {
				hash = hash << 2;
				hash ^= CGAT_REVERSE_JUMP_TABLE[buffer.get(k - i - 1)];
			}
		} else {
			for (int i = 0; i < k; i++) {
				hash = hash << 2;
				hash ^= CGAT_JUMP_TABLE[buffer.get(i)];
			}
		}
		return hash;
	}

	private void putViaHash(long hash) {
		int nextHash;
		int index;
		int[][] perm;		
		
		for (int i = 0; i < t.length; i++) {
			perm = t[i];
			nextHash = 0;
			for (int j = 0; j < perm.length; j++) {
				nextHash = nextHash + perm[j][(int) (hash & 0b11111111)];
				hash >>>= 8;
			}
			nextHash = nextHash % sizeInBits;
			if (nextHash < 0) {
				nextHash = -nextHash;
			}
						
			index = nextHash >>> 6;
			bits[index] = bits[index] | (1L << (nextHash & 0b111111));
		}
	}
	
	public boolean contains(byte[] seq, int start, boolean reverse, int[] badPos) {
		long hash = hash(seq, start, reverse);
		return containsHash(hash);
	}

	private boolean containsHash(long hash) {
		int nextHash;
		int index;
		int[][] perm;		
		
		for (int i = 0; i < t.length; i++) {
			perm = t[i];
			nextHash = 0;
			for (int j = 0; j < perm.length; j++) {
				nextHash = nextHash + perm[j][(int) (hash & 0b11111111)];
				hash >>>= 8;
			}
			nextHash = nextHash % sizeInBits;
			if (nextHash < 0) {
				nextHash = -nextHash;
			}
						
			index = nextHash >>> 6;
			if (((bits[index] >> (nextHash & 0b111111)) & 1L) == 0) {
				return false;
			}
		}
		return true;
	}

	public boolean contains(CGATRingBuffer buffer, boolean reverse, int[] badPos) {
		long hash = hash(buffer, reverse);
		return containsHash(hash);
	}
}