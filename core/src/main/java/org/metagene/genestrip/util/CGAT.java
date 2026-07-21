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

/**
 * Low-level helpers for encoding DNA bases (C, G, A, T) into 2-bit codes and k-mers into longs,
 * as well as decoding, complementing and reverse-complementing them.
 */
public class CGAT {
	// Static utility class - not meant to be instantiated.
	private CGAT() {
	}

	/** Lookup table mapping each nucleotide byte to its Watson-Crick complement; non-CGAT entries are -1. */
	private static final byte[] CGAT_COMPLEMENT = new byte[127];
	/** Lookup table mapping each nucleotide byte to its straight 2-bit code (C=0, G=1, A=2, T=3); non-CGAT entries are -1. */
	static final int[] CGAT_JUMP_TABLE;
	/** Lookup table mapping each nucleotide byte to the 2-bit code of its complement (C=1, G=0, A=3, T=2); non-CGAT entries are -1. */
	static final int[] CGAT_REVERSE_JUMP_TABLE;

	/** Maps a 2-bit code (0-3) back to its nucleotide letter C, G, A or T. */
	private static final byte[] DECODE_TABLE = new byte[] { 'C', 'G', 'A', 'T' };

	/** Per-k bit masks retaining the low 2k bits, used when appending a base to a straight k-mer encoding. */
	static final long SHIFT_FILTERS_STRAIGHT[] = new long[33]; // TODO ~(-1 << ((k - 1) * 2))
	/** Per-k shift amounts ((k-1)*2), used when prepending a base to a reverse k-mer encoding. */
	static final long SHIFT_FILTERS_REVERSE[] = new long[33]; // TODO ((k - 1) * 2)

	static {
		for (int i = 0; i < CGAT_COMPLEMENT.length; i++) {
			CGAT_COMPLEMENT[i] = -1;
		}
		CGAT_COMPLEMENT['C'] = 'G';
		CGAT_COMPLEMENT['G'] = 'C';
		CGAT_COMPLEMENT['A'] = 'T';
		CGAT_COMPLEMENT['T'] = 'A';

		CGAT_JUMP_TABLE = new int[Byte.MAX_VALUE + 1];
		CGAT_REVERSE_JUMP_TABLE = new int[Byte.MAX_VALUE + 1];
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

	/**
	 * Returns a fresh copy of the 2-bit-code-to-nucleotide decode table (index 0-3 maps to C, G, A, T).
	 * A copy is handed out so callers cannot mutate the shared internal table; fetch it once and reuse it.
	 *
	 * @return a new {@code byte[]} copy of the decode table
	 */
	public static byte[] newDecodeTable() {
		return DECODE_TABLE.clone();
	}

	/**
	 * Upper-cases only the four nucleotide letters; every other byte is returned unchanged.
	 *
	 * @param c the byte to upper-case
	 * @return the uppercase nucleotide letter, or {@code c} unchanged if it is not a, c, g or t
	 */
	public static byte cgatToUpperCase(byte c) {
		switch (c) {
		case 'a': return 'A';
		case 'c': return 'C';
		case 'g': return 'G';
		case 't': return 'T';
		default:  return c;
		}
	}

	/**
	 * Returns whether the given byte is one of the uppercase nucleotide letters C, G, A or T.
	 *
	 * @param c the byte to test
	 * @return {@code true} if {@code c} is C, G, A or T
	 */
	public static boolean isCGAT(byte c) {
		return c == 'C' || c == 'G' || c == 'A' || c == 'T';
	}

	/**
	 * Returns the Watson-Crick complement of the given nucleotide byte, or -1 for non-CGAT input.
	 *
	 * @param c the nucleotide byte to complement
	 * @return the complement byte, or -1 if {@code c} is not a CGAT base
	 */
	public static byte toComplement(byte c) {
		return CGAT_COMPLEMENT[c];
	}

	/**
	 * Encodes the k bases starting at {@code start} into the canonical (standard) k-mer, i.e. the
	 * larger of the straight and reverse-complement 2-bit encodings.
	 *
	 * @param seq    the byte array holding the bases
	 * @param start  the index of the first base to encode
	 * @param k      the number of bases (k-mer length)
	 * @param badPos if non-null, {@code badPos[0]} is set to the position of the first non-CGAT base,
	 *               or -1 if none.
	 * @return the canonical k-mer encoding, or -1 if a non-CGAT base is encountered.
	 */
	public static long kMerToLong(byte[] seq, int start, int k, int[] badPos) {
		long reverse = kMerToLongReverse(seq, start, k, badPos);
		long straight = kMerToLongStraight(seq, start, k, badPos);
		return standardKMer(reverse, straight);
	}

	/**
	 * Returns the canonical k-mer, i.e. the larger of the two given encodings.
	 *
	 * @param straight the straight k-mer encoding
	 * @param reverse  the reverse-complement k-mer encoding
	 * @return the larger of {@code straight} and {@code reverse}
	 */
	public static long standardKMer(long straight, long reverse) {
		return straight > reverse ? straight : reverse;
	}

	/**
	 * Encodes the k bases starting at {@code start} into a 2-bit-per-base long in reading direction.
	 *
	 * @param seq    the byte array holding the bases
	 * @param start  the index of the first base to encode
	 * @param k      the number of bases (k-mer length)
	 * @param badPos if non-null, {@code badPos[0]} is set to the position of the first non-CGAT base,
	 *               or -1 if none.
	 * @return the straight k-mer encoding, or -1 if a non-CGAT base is encountered.
	 */
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
	
	/**
	 * Decodes a straight 2-bit-per-base k-mer encoding back into k nucleotide bytes, written to
	 * {@code res} starting at {@code start}.
	 *
	 * @param kmer  the straight k-mer encoding to decode
	 * @param res   the byte array to write the decoded bases into
	 * @param start the index in {@code res} at which to write the first base
	 * @param k     the number of bases (k-mer length)
	 */
	public static void longToKMerStraight(long kmer, byte[] res, int start, int k) {
		for (int i = k - 1; i >= 0; i--) {
			// & 3 and >>> 2 (not % 4 and >> 2) so a k-mer whose top bit is set decodes correctly.
			res[start + i] = DECODE_TABLE[(int) (kmer & 3)];
			kmer = kmer >>> 2;
		}
	}

	/**
	 * Returns the straight k-mer encoding obtained by appending base {@code bp} to the given k-mer
	 * and dropping its oldest base, or -1 if {@code bp} is not a CGAT base.
	 *
	 * @param kmer the current straight k-mer encoding
	 * @param bp   the nucleotide byte to append
	 * @param k    the number of bases (k-mer length)
	 * @return the updated straight k-mer encoding, or -1 if {@code bp} is not a CGAT base
	 */
	public static long nextKMerStraight(final long kmer, final byte bp, final int k) {
		int c = CGAT_JUMP_TABLE[bp];
		if (c == -1) {
			return -1L;
		}
		return ((kmer << 2) & SHIFT_FILTERS_STRAIGHT[k]) | (long) c;
	}

	/**
	 * Returns the reverse-complement k-mer encoding obtained by prepending the complement of base
	 * {@code bp} to the given reverse k-mer and dropping its oldest base, or -1 if {@code bp} is not
	 * a CGAT base.
	 *
	 * @param kmer the current reverse-complement k-mer encoding
	 * @param bp   the nucleotide byte whose complement to prepend
	 * @param k    the number of bases (k-mer length)
	 * @return the updated reverse-complement k-mer encoding, or -1 if {@code bp} is not a CGAT base
	 */
	public static long nextKMerReverse(final long kmer, final byte bp, final int k) {
		int c = CGAT_REVERSE_JUMP_TABLE[bp];
		if (c == -1) {
			return -1L;
		}
		return (kmer >>> 2) | (((long) c) << SHIFT_FILTERS_REVERSE[k]);
	}

	/**
	 * Encodes the reverse complement of the k bases starting at {@code start} into a 2-bit-per-base
	 * long.
	 *
	 * @param seq    the byte array holding the bases
	 * @param start  the index of the first base to encode
	 * @param k      the number of bases (k-mer length)
	 * @param badPos if non-null, {@code badPos[0]} is set to the position of the first non-CGAT base,
	 *               or -1 if none.
	 * @return the reverse-complement k-mer encoding, or -1 if a non-CGAT base is encountered.
	 */
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

	/**
	 * Reverse-complements the entire array in place.
	 *
	 * @param seq the byte array to reverse-complement
	 */
	public static void reverse(byte[] seq) {
		reverse(seq, 0, seq.length);
	}

	/**
	 * Reverse-complements seq[start, start+k) in place (k is the region length, not an end index).
	 *
	 * @param seq   the byte array to reverse-complement
	 * @param start the index of the first base of the region
	 * @param k     the length of the region to reverse-complement
	 */
	public static void reverse(byte[] seq, int start, int k) {
		int end = start + k - 1;
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
