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

public class CGATLongBuffer {
	protected final int size;
	
	protected int bpCounter;
	private long kmer;
	private long reverseKmer;
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
	protected final int srl0Buffer[];
	protected final int srl1Buffer[];
	protected final int srl2Buffer[];

	public CGATLongBuffer(int size) {
		this(size, -1);
	}

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
	
	public final long getKMer() {
		return filled ? kmer : -1;
	}

	public final long getReverseKMer() {
		return filled ? reverseKmer : -1;
	}

	public final long getStandardKMer() {
		if (filled) {
			return CGAT.standardKMer(kmer, reverseKmer);
		}
		return -1;
	}

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

	public final boolean isFilled() {
		return filled;
	}

	public final int getSize() {
		return size;
	}

	public final boolean isDust() {
		return d > maxDust;
	}

	public final int getDustValue() {
		return d;
	}
	
	public final int getMaxDust() {
		return maxDust;
	}
}