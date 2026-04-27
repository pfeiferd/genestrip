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
package org.metagene.genestrip.util;

public class CGAT {
	public static final byte[] CGAT_TO_UPPER_CASE = new byte[256];
	public static final byte[] CGAT_COMPLEMENT = new byte[127];
	public static final int[] CGAT_JUMP_TABLE;
	public static final int[] CGAT_REVERSE_JUMP_TABLE;

	public static final byte[] DECODE_TABLE = new byte[] { 'C', 'G', 'A', 'T' };
	public static final byte[] REVERSE_DECODE_TABLE = new byte[] { 'G', 'C', 'T', 'A' };

	public static long SHIFT_FILTERS_STRAIGHT[] = new long[33]; // TODO ~(-1 << ((k - 1) * 2))
	public static long SHIFT_FILTERS_REVERSE[] = new long[33]; // TODO ((k - 1) * 2)

	static {
		for (int i = 0; i < CGAT_TO_UPPER_CASE.length; i++) {
			CGAT_TO_UPPER_CASE[i] = (byte) (i - 128);
		}
		CGAT_TO_UPPER_CASE[128 + 'c'] = 'C';
		CGAT_TO_UPPER_CASE[128 + 'g'] = 'G';
		CGAT_TO_UPPER_CASE[128 + 'a'] = 'A';
		CGAT_TO_UPPER_CASE[128 + 't'] = 'T';

		for (int i = 0; i < CGAT_COMPLEMENT.length; i++) {
			CGAT_COMPLEMENT[i] = -1;
		}
		CGAT_COMPLEMENT['C'] = 'G';
		CGAT_COMPLEMENT['G'] = 'C';
		CGAT_COMPLEMENT['A'] = 'T';
		CGAT_COMPLEMENT['T'] = 'A';

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

		for (int i = 0; i < SHIFT_FILTERS_STRAIGHT.length - 1; i++) {
			SHIFT_FILTERS_STRAIGHT[i] = ~(-1L << (i * 2));
		}
		SHIFT_FILTERS_STRAIGHT[32] = -1;
		for (int i = 1; i < SHIFT_FILTERS_REVERSE.length; i++) {
			SHIFT_FILTERS_REVERSE[i] = (i - 1) * 2;
		}
	}

	public static byte cgatToUpperCase(byte c) {
		return CGAT_TO_UPPER_CASE[128 + c];
	}

	public static boolean isCGAT(byte c) {
		return c == 'C' || c == 'G' || c == 'A' || c == 'T';
	}

	public static byte toComplement(byte c) {
		return CGAT_COMPLEMENT[c];
	}

	public static long kMerToLong(byte[] seq, int start, int k, int[] badPos) {
		long reverse = kMerToLongReverse(seq, start, k, badPos);
		long straight = kMerToLongStraight(seq, start, k, badPos);
		return standardKMer(reverse, straight);
	}

	public static long standardKMer(long straight, long reverse) {
		return straight > reverse ? straight : reverse;
	}

	public static long kMerToLongStraight(final byte[] seq, final int start, final int k, final int[] badPos) {
		long res = 0;
		int c;
		if (badPos != null) {
			badPos[0] = -1;
		}
		int max = start + k;
		for (int i = start; i < max; i++) {
			// Inlined: res = Long.rotateLeft(res, 2);
			res = (res << 2) | (res >>> -2);
			c = CGAT_JUMP_TABLE[seq[i]];
			if (c == -1) {
				if (badPos != null) {
					badPos[0] = i;
				}
				return -1L;
			}
			res += c;
		}

		return res;
	}
	
	public static void longToKMerStraight(long kmer, byte[] res, int start, int k) {
		assert(k <= 31);
		for (int i = k - 1; i >= 0; i--) {
			res[start + i] = DECODE_TABLE[(int) (kmer % 4)]; 
			kmer = kmer >> 2;
		}
	}

	public static long nextKMerStraight(final long kmer, final byte bp, final int k) {
		int c = CGAT_JUMP_TABLE[bp];
		if (c == -1) {
			return -1L;
		}
		return ((kmer << 2) & SHIFT_FILTERS_STRAIGHT[k]) | (long) c;
	}

	public static long nextKMerReverse(final long kmer, final byte bp, final int k) {
		int c = CGAT_REVERSE_JUMP_TABLE[bp];
		if (c == -1) {
			return -1L;
		}
		return (kmer >>> 2) | (((long) c) << SHIFT_FILTERS_REVERSE[k]);
	}

	public static long kMerToLongReverse(final byte[] seq, final int start, final int k, final int[] badPos) {
		long res = 0;
		int c;
		if (badPos != null) {
			badPos[0] = -1;
		}
		for (int i = start + k - 1; i >= start; i--) {
			// Inlined: res = Long.rotateLeft(res, 2);
			res = (res << 2) | (res >>> -2);
			c = CGAT_REVERSE_JUMP_TABLE[seq[i]];
			if (c == -1) {
				if (badPos != null) {
					badPos[0] = i;
				}
				return -1L;
			}
			res += c;
		}

		return res;
	}

	public static void reverse(byte[] seq) {
		reverse(seq, 0, seq.length);
	}

	public static void reverse(byte[] seq, int start, int k) {
		int end = k - 1;
		byte h;
		while (start < end) {
			h = CGAT_COMPLEMENT[seq[start]];
			seq[start] = CGAT_COMPLEMENT[seq[end]];
			seq[end] = h;
			start++;
			end--;
		}
		if (start == end) {
			seq[start] = CGAT_COMPLEMENT[seq[start]];
		}
	}
}
