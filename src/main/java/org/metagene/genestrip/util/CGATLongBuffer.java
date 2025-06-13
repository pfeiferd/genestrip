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
	private final int seqMarks[];
	protected final int size;
	
	protected int bpCounter;
	private long kmer;
	private long reverseKmer;
	protected boolean filled;

	private final int maxDust;
	private final int[] dustFunction;

	private byte lastChar;
	private int sumDust;
	private int seqCount;

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
			dustFunction = new int[size];
			initDustFunction(dustFunction);
			seqMarks = new int[size];
		}
		else {
			seqMarks = null;
			dustFunction = null;
		}
		reset();
	}
	
	protected void initDustFunction(int[] dustFunction) {		
		if (dustFunction.length > 0) {
			dustFunction[0] = 1;
		}
		if (dustFunction.length > 1) {
			dustFunction[1] = 1;
		}
		if (dustFunction.length > 2) {
			dustFunction[2] = 1;
		}
		for (int i = 3; i < dustFunction.length; i++) {
			dustFunction[i] = dustFunction[i - 1] + dustFunction[i - 2];
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
				if (c == lastChar) {
					seqMarks[(bpCounter - seqCount + size) % size]++;
					sumDust += dustFunction[seqCount];
					if (seqCount < size - 1) {
						seqCount++;						
					}
				} else {
					seqCount = 0;
				}
				lastChar = c;
			}
			bpCounter++;
			if (bpCounter == size) {
				bpCounter = 0;
				filled = true;
			}
			if (filled && maxDust >= 0) {
				int oldCount = seqMarks[bpCounter];
				seqMarks[bpCounter] = 0;
				if (oldCount > 0) {
					sumDust -= dustFunction[oldCount - 1];
					seqMarks[(bpCounter + 1) % size] = oldCount - 1;
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
		seqCount = 0;
		sumDust = maxDust >= 0 ? 0 : -1;
		
		if (seqMarks != null) {
			for (int i = 0; i < size; i++) {
				seqMarks[i] = 0;
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
		return sumDust > maxDust;
	}

	public final int getDustValue() {
		return sumDust;
	}
	
	public final int getMaxDust() {
		return maxDust;
	}
}