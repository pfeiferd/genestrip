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

import org.apache.commons.codec.digest.MurmurHash3;

import it.unimi.dsi.fastutil.BigArrays;

public class LargeMurmurCGATBloomFilter extends TwoLongsCGATBloomFilter implements Serializable {
	private static final long serialVersionUID = 1L;
	
	private final int[] hashFactors1;
	private final int[] hashFactors2;
	protected long[][] bits;
	
	public LargeMurmurCGATBloomFilter(int k, long expectedInsertions, double fpp) {
		super(k, expectedInsertions, fpp);
		
		Random r = new Random(42);
		hashFactors1 = new int[hashes];
		hashFactors2 = new int[hashes];
		for (int i = 0; i < hashFactors1.length; i++) {
			hashFactors1[i] = r.nextInt();
			hashFactors2[i] = r.nextInt();
		}
	}
		
	@Override
	protected void clearArray() {
		BigArrays.fill(bits, 0);
	}
		
	@Override
	protected void initBitArray() {
		if (bits == null) {
			bits = BigArrays.ensureCapacity(BigArrays.wrap(new long[0]), size);
		}
		else {
			bits = BigArrays.ensureCapacity(bits, size);
		}
	}

	@Override
	protected boolean containsHash(long hash1, long hash2) {
		long hash;
		long index;

		for (int i = 0; i < hashes; i++) {
			hash = combineLongHashesLarge(hash1, hash2, i);
			index = ((hash >>> 6) % size);
			if (((BigArrays.get(bits, index) >> (hash & 0b111111)) & 1L) == 0) {
				return false;
			}
		}
		return true;
	}
	
	@Override
	protected void putViaHash(long hash1, long hash2) {
		long hash;
		long index;

		for (int i = 0; i < hashes; i++) {
			hash = combineLongHashesLarge(hash1, hash2, i);
			index = ((hash >>> 6) % size);
			BigArrays.set(bits, index, BigArrays.get(bits, index) | (1L << (hash & 0b111111)));
		}
	}
	
	protected long combineLongHashesLarge(long hash1, long hash2, int i) {
		long a = MurmurHash3.hash32(hash1, hash2, hashFactors1[i]);
		int b = MurmurHash3.hash32(hash1, hash2, hashFactors2[i]);
		
		return a << 32 | b;
	}
	
	@Override
	protected int combineLongHashes(long hash1, long hash2, int i) {
		throw new UnsupportedOperationException("Not to be called in this variant.");
	}
}