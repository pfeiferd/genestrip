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

import static org.metagene.genestrip.util.CGAT.*;

/**
 * A streaming ring buffer that incrementally computes the straight and reverse-complement 2-bit
 * encodings of the last k bases fed to it via {@link #put(byte)}, optionally tracking a DUST
 * low-complexity score so that low-complexity k-mers can be filtered out.
 */
public class CGATLongBuffer {
	/** The k-mer length in bases. */
	protected final int size;

	/** The current write position within the ring buffer. */
	protected int bpCounter;
	private long kmer;
	private long reverseKmer;
	/** Whether at least {@code size} bases have been buffered. */
	protected boolean filled;

	private final int maxDust;
	private final int[] dustFunctionDiff;

	private byte l1Char;
	private byte l2Char;
	private byte l3Char;
	private int d;
	private int srl0;
	private int srl1;
	private int srl2;
	/** Ring buffer of run lengths for the DUST computation at offset 0. */
	protected final int srl0Buffer[];
	/** Ring buffer of run lengths for the DUST computation at offset 1. */
	protected final int srl1Buffer[];
	/** Ring buffer of run lengths for the DUST computation at offset 2. */
	protected final int srl2Buffer[];

	/**
	 * Creates a buffer for k-mers of the given length with DUST filtering disabled.
	 *
	 * @param size the k-mer length in bases (at most 32)
	 */
	public CGATLongBuffer(int size) {
		this(size, -1);
	}

	/**
	 * Creates a buffer for k-mers of the given length (at most 32).
	 *
	 * @param size    the k-mer length in bases (at most 32)
	 * @param maxDust the maximum tolerated DUST score, or a negative value to disable DUST tracking.
	 */
	public CGATLongBuffer(int size, int maxDust) {
		if (size > 32) {
			throw new IllegalArgumentException("size must be <= 32");
		}
		if (maxDust > Short.MAX_VALUE) {
			throw new IllegalArgumentException("Unreasonable size for maxDust :" + maxDust);
		}
		this.size = size;
		this.maxDust = maxDust;
		if (maxDust >= 0) {
			dustFunctionDiff = new int[size];
			initDustFunctionDiff(dustFunctionDiff);
			srl0Buffer = new int[size];
			srl1Buffer = new int[size];
			srl2Buffer = new int[size];
		}
		else {
			srl0Buffer = srl1Buffer = srl2Buffer = null;
			dustFunctionDiff = null;
		}
		reset();
	}
	
	/**
	 * Initializes the DUST weighting table used to incrementally update the DUST score.
	 *
	 * @param dustFunctionDiff the array to fill with the DUST weighting differences
	 */
	protected void initDustFunctionDiff(int[] dustFunctionDiff) {
		// For the streaming approach we must work with
		// the differences fib(n + 1) - fib(n - 1) - which are slightly different from the
		// originally used Fibonacci function with fib(0) = 0, fib(1) = 1, fib(2) = 2, fib(n) = fib(n - 1) + fib(n - 2)
		if (dustFunctionDiff.length > 0) {
			dustFunctionDiff[0] = 1;
		}
		if (dustFunctionDiff.length > 1) {
			dustFunctionDiff[1] = 1;
		}
		if (dustFunctionDiff.length > 2) {
			dustFunctionDiff[2] = 1;
		}
		for (int i = 3; i < dustFunctionDiff.length; i++) {
			dustFunctionDiff[i] = dustFunctionDiff[i - 1] + dustFunctionDiff[i - 2];
		}
	}

	/**
	 * Feeds the next base into the buffer, updating the straight and reverse-complement encodings and
	 * (if enabled) the DUST score.
	 *
	 * @param c the next base to feed into the buffer
	 * @return the updated straight k-mer encoding, or -1 if {@code c} is not a CGAT base, in which
	 *         case the buffer is reset.
	 */
	public final long put(byte c) {
		// This is inlined: kmer = CGAT.nextKMerStraight(kmer, c, size);
		// And also: reverseKmer = CGAT.nextKMerReverse(reverseKmer, c, size)
		int bp = CGAT_JUMP_TABLE[c]; // Inlined.
		if (bp == -1) {  // Inlined.
			reset();
			return -1L;
		}
		else {
			kmer = ((kmer << 2) & SHIFT_FILTERS_STRAIGHT[size]) | (long) bp;  // Inlined.
			reverseKmer = (reverseKmer >>> 2) | (((long) CGAT_REVERSE_JUMP_TABLE[c]) << SHIFT_FILTERS_REVERSE[size]);
			if (maxDust >= 0) {
				if (c == l1Char) {
					int pos = bpCounter - 1 - srl0;
					if (pos < 0) {
						pos += size;
					}
					srl0Buffer[pos]++;
					d += dustFunctionDiff[srl0];
					if (srl0 < size - 1) {
						srl0++;
					}
				} else {
					srl0 = 0;
				}
				if (c == l2Char) {
					int pos = bpCounter - 2 - srl1;
					if (pos < 0) {
						pos += size;
					}
					srl1Buffer[pos]++;
					d += dustFunctionDiff[srl1];
					if (srl1 < size - 2) {
						srl1++;
					}
				} else {
					srl1 = 0;
				}
				if (c == l3Char) {
					int pos = bpCounter - 3 - srl2;
					if (pos < 0) {
						pos += size;
					}
					srl2Buffer[pos]++;
					d += dustFunctionDiff[srl2];
					if (srl2 < size - 3) {
						srl2++;
					}
				} else {
					srl2 = 0;
				}
				l3Char = l2Char;
				l2Char = l1Char;
				l1Char = c;
			}
			int oldBp = bpCounter;
			bpCounter++;
			if (bpCounter == size) {
				bpCounter = 0;
				filled = true;
			}
			if (filled && maxDust >= 0) {
				int oldCount = srl0Buffer[oldBp];
				srl0Buffer[oldBp] = 0;
				if (oldCount > 0) {
					d -= dustFunctionDiff[oldCount - 1];
					srl0Buffer[bpCounter] = oldCount - 1;
				}

				oldCount = srl1Buffer[oldBp];
				srl1Buffer[oldBp] = 0;
				if (oldCount > 0) {
					d -= dustFunctionDiff[oldCount - 1];
					srl1Buffer[bpCounter] = oldCount - 1;
				}

				oldCount = srl2Buffer[oldBp];
				srl2Buffer[oldBp] = 0;
				if (oldCount > 0) {
					d -= dustFunctionDiff[oldCount - 1];
					srl2Buffer[bpCounter] = oldCount - 1;
				}
			}
			return kmer;
		}
	}
	
	/**
	 * Returns the current straight k-mer encoding, or -1 if fewer than k bases have been buffered.
	 *
	 * @return the straight k-mer encoding, or -1 if not yet filled
	 */
	public final long getKMer() {
		return filled ? kmer : -1;
	}

	/**
	 * Returns the current reverse-complement k-mer encoding, or -1 if fewer than k bases have been
	 * buffered.
	 *
	 * @return the reverse-complement k-mer encoding, or -1 if not yet filled
	 */
	public final long getReverseKMer() {
		return filled ? reverseKmer : -1;
	}

	/**
	 * Returns the current canonical (standard) k-mer encoding, or -1 if fewer than k bases have been
	 * buffered.
	 *
	 * @return the canonical k-mer encoding, or -1 if not yet filled
	 */
	public final long getStandardKMer() {
		if (filled) {
			return CGAT.standardKMer(kmer, reverseKmer);
		}
		return -1;
	}

	/**
	 * Clears all buffered bases and DUST state so that the buffer starts empty again.
	 */
	public final void reset() {
		bpCounter = 0;
		kmer = 0;
		reverseKmer = 0;
		filled = false;
		srl0 = 0;
		srl1 = 0;
		srl2 = 0;
		l1Char = -1;
		l2Char = -1;
		l3Char = -1;
		d = maxDust >= 0 ? 0 : -1;
		
		if (maxDust >= 0) {
			for (int i = 0; i < size; i++) {
				srl0Buffer[i] = 0;
				srl1Buffer[i] = 0;
				srl2Buffer[i] = 0;
			}
		}
	}

	/**
	 * Returns whether at least {@code size} bases have been buffered.
	 *
	 * @return {@code true} if the buffer is filled
	 */
	public final boolean isFilled() {
		return filled;
	}

	/**
	 * Returns the k-mer length in bases.
	 *
	 * @return the k-mer length
	 */
	public final int getSize() {
		return size;
	}

	/**
	 * Returns whether the current k-mer's DUST score exceeds the configured maximum (i.e. it is
	 * low-complexity).
	 *
	 * @return {@code true} if the current k-mer is low-complexity
	 */
	public final boolean isDust() {
		return d > maxDust;
	}

	/**
	 * Returns the current DUST score, or -1 if DUST tracking is disabled.
	 *
	 * @return the current DUST score
	 */
	public final int getDustValue() {
		return d;
	}

	/**
	 * Returns the maximum tolerated DUST score, or a negative value if DUST tracking is disabled.
	 *
	 * @return the maximum tolerated DUST score
	 */
	public final int getMaxDust() {
		return maxDust;
	}
}