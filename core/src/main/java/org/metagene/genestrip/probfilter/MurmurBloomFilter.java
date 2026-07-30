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
package org.metagene.genestrip.probfilter;

/**
 * A {@link BloomFilter} that hashes k-mers with MurmurHash3 rather than with the bit-mixing
 * finalizer of {@link BloomFilter#hash(long, int)}. Both mix, so this inherits the multiply-shift
 * {@link BloomFilter#reduce(long)} unchanged; only the hash differs. The {@code XOR} variant of the
 * family, which does not mix and hence reduces by a modulo, is {@link XORBloomFilter}.
 */
public class MurmurBloomFilter extends BloomFilter {
	// Bumped from 1: the inherited reduce() changed from a modulo to a multiply-shift, so a key maps to
	// a different bit than before. A filter serialized by an older version would answer queries wrongly
	// rather than merely differently, hence it must fail to load instead - such filters have to be
	// regenerated. XORBloomFilter keeps its modulo and hence its serial version.
	private static final long serialVersionUID = 2L;

	/**
	 * Creates a filter with the given target false-positive probability, sized for
	 * {@code expectedInsertions} k-mers.
	 *
	 * @param fpp the target false-positive probability
	 * @param expectedInsertions the expected number of k-mers to be inserted
	 */
	public MurmurBloomFilter(double fpp, long expectedInsertions) {
		super(fpp, expectedInsertions);
	}

	@Override
	protected final long hash(final long data, final int i) {
		return MurmurHash3DropIn.hash64(data, hashFactors[i]);
	}

	/**
	 * A drop-in reimplementation of the 64-bit MurmurHash3 hash (as in Apache Commons Codec) that,
	 * unlike the deprecated original, supports a custom seed.
	 */
	public static class MurmurHash3DropIn {
		// Static utility class - not meant to be instantiated.
		private MurmurHash3DropIn() {
		}

		/*
		 * The following code was copied from
		 * org.apache.commons.codec.digest.MurmurHash3.hash64()
		 *
		 * The original function is deprecated and does not support changing seeds.
		 * That's why its implemented here...
		 */

		// Constants for 128-bit variant
		private static final long C1 = 0x87c37b91114253d5L;
		private static final long C2 = 0x4cf5ad432745937fL;
		private static final int R1 = 31;
		private static final int R2 = 27;
		private static final int M = 5;
		private static final int N1 = 0x52dce729;

		/**
		 * Computes the 64-bit MurmurHash3 hash of the given 8-byte value, using {@code hashBase} as the
		 * seed.
		 *
		 * @param data     the 8-byte value to hash
		 * @param hashBase the seed value
		 * @return the 64-bit MurmurHash3 hash
		 */
		public static final long hash64(final long data, long hashBase) {
			long hash = hashBase;

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
}