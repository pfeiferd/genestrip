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

public abstract class TwoLongsCGATBloomFilter extends AbstractCGATBloomFilter implements Serializable {
	private static final long serialVersionUID = 1L;

	public TwoLongsCGATBloomFilter(int k, long expectedInsertions, double fpp) {
		super(k, expectedInsertions, fpp);
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

	private long hash(byte[] seq, int start, boolean even, boolean reverse) {
		long hash = 0;
		int c;
		
		if (reverse) {
			for (int i = even ? 0 : 1; i < k; i += 2) {
				hash = hash << 2;
				c = CGAT_REVERSE_JUMP_TABLE[seq[k - i - 1]];
				if (c == -1) {
					return 0;
				}
				hash += c;
			}
		} else {
			for (int i = even ? 0 : 1; i < k; i += 2) {
				hash = hash << 2;
				c = CGAT_JUMP_TABLE[seq[i]];
				if (c == -1) {
					return 0;
				}
				hash += c;				
			}
		}
		return hash;
	}

	private long hash(CGATRingBuffer buffer, boolean even, boolean reverse) {
		assert (buffer.getSize() == k);

		long hash = 0;
		int c;
		
		if (reverse) {
			for (int i = even ? 0 : 1; i < k; i += 2) {
				hash = hash << 2;
				c = CGAT_REVERSE_JUMP_TABLE[buffer.get(k - i - 1)];
				if (c == -1) {
					return 0;
				}
				hash += c;				
			}
		} else {
			for (int i = even ? 0 : 1; i < k; i += 2) {
				hash = hash << 2;
				c = CGAT_JUMP_TABLE[buffer.get(i)];
				if (c == -1) {
					return 0;
				}
				hash += c;				
			}
		}
		return hash == 0 ? 1 : hash;
	}

	private void putViaHash(long hash1, long hash2) {
		long hash;
		int index;
		
		for (int i = 0; i < hashes; i++) {
			hash = combineLongHashes(hash1, hash2, i);
			index = (int) ((hash >>> 6) % bits.length);
			bits[index] = bits[index] | (1L << (hash & 0b111111L));
		}
	}
	
	protected abstract long combineLongHashes(long hash1, long hash2, int i);

	public boolean contains(byte[] seq, int start, boolean reverse) {
		long hash1 = hash(seq, start, true, reverse);
		if (hash1 == 0) {
			return false;
		}
		long hash2 = hash(seq, start, false, reverse);
		if (hash2 == 0) {
			return false;
		}
		return containsHash(hash1, hash2);
	}

	private boolean containsHash(long hash1, long hash2) {
		long hash;
		int index;
		
		for (int i = 0; i < hashes; i++) {
			hash = combineLongHashes(hash1, hash2, i);
			index = (int) ((hash >>> 6) % bits.length);
			if (((bits[index] >> (hash & 0b111111L)) & 1L) == 0) {
				return false;
			}
		}
		return true;
	}

	public boolean contains(CGATRingBuffer buffer, boolean reverse) {
		long hash1 = hash(buffer, true, reverse);
		if (hash1 == 0) {
			return false;
		}
		long hash2 = hash(buffer, false, reverse);
		if (hash2 == 0) {
			return false;
		}
		return containsHash(hash1, hash2);
	}
}