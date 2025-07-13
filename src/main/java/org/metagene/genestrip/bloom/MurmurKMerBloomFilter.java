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

public class MurmurKMerBloomFilter extends AbstractKMerBloomFilter {
	public MurmurKMerBloomFilter(int k, double fpp) {
		super(k, fpp);
	}

	/*
	public final boolean containsViaHashInlined(final long data) {
		for (int i = 0; i < hashes; i++) {
			long hash = hashFactors[i];

			// Reverse byte inlined from Long.reverseBytes()
			long k1 = (data & 0x00ff00ff00ff00ffL) << 8 | (data >>> 8) & 0x00ff00ff00ff00ffL;
			k1 = (k1 << 48) | ((k1 & 0xffff0000L) << 16) | ((k1 >>> 16) & 0xffff0000L) | (k1 >>> 48);

			final int length = Long.BYTES;
			// mix functions
			k1 *= 0x87c37b91114253d5L;
			// Rotate left inlined from Long.rotateLeft()
			k1 = (k1 << 31) | (k1 >>> -31);
			k1 *= C2;
			hash ^= k1;
			// Rotate left inlined from Long.rotateLeft()
			hash = ((hash << 27) | (hash >>> -27)) * 5 + 0x52dce729;
			// finalization
			hash ^= length;

			// Inlined from MurmurHash3.fmix64()
			hash ^= (hash >>> 33);
			hash *= 0xff51afd7ed558ccdL;
			hash ^= (hash >>> 33);
			hash *= 0xc4ceb9fe1a85ec53L;
			hash ^= (hash >>> 33);
			if (bitVector.largeBits != null) {
				long arrayIndex = ((hash >>> 6) % bitVector.size);
				if (((bitVector.largeBits[(int) (arrayIndex >>> 27)][(int) (arrayIndex & 134217727)] >> (hash & 0b111111)) & 1L) != 1) {
					return false;
				}
			} else {
				int arrayIndex = (int) ((hash >>> 6) % bitVector.size);
				if ((((bitVector.bits[arrayIndex] >> (hash & 0b111111)) & 1L) != 1)) {
					return false;
				}
			}
		}
		return true;
	}
	 */

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

	protected final long hash(final long data, final int i) {
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

		return hash ^ data;
	}
}