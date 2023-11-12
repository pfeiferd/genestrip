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

import java.io.Serializable;

public class CGATLongBuffer implements Serializable {
	private static final long serialVersionUID = 1L;

	private final int seqMarks[];
	private final int size;
	
	private int bpCounter;
	private long kmer;
	private boolean filled;

	private final int maxDust;
	private final int[] fib;

	private byte lastChar;
	private int sumDust;
	private int seqCount;

	public CGATLongBuffer(int size) {
		this(size, -1);
	}

	public CGATLongBuffer(int size, int maxDust) {
		if (maxDust > Short.MAX_VALUE) {
			throw new IllegalArgumentException("Unreasonable size for maxDust :" + maxDust);
		}
		this.size = size;
		this.maxDust = maxDust;
		if (maxDust >= 0) {
			fib = new int[size];
			if (size > 0) {
				fib[0] = 1;
			}
			if (size > 1) {
				fib[1] = 1;
			}
			if (size > 2) {
				fib[2] = 1;
			}
			for (int i = 3; i < fib.length; i++) {
				fib[i] = fib[i - 1] + fib[i - 2];
			}
			seqMarks = new int[size];
		}
		else {
			seqMarks = null;
			fib = null;
		}
		reset();
	}

	public void put(byte c) {
		kmer = CGAT.nextKMerStraight(kmer, c, size);
		if (kmer == -1) {
			reset();
		}
		else {
			if (maxDust >= 0) {
				if (c == lastChar) {
					seqMarks[bpCounter - seqCount]++;
					sumDust += fib[seqCount];
					if (seqCount < size) {
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
					sumDust -= fib[oldCount - 1];
					seqMarks[(bpCounter + 1) % size] = oldCount - 1;
				}
			}
		}
	}
	
	public long getKMer() {
		return filled ? kmer : -1;
	}

	public void reset() {
		bpCounter = 0;
		kmer = 0;
		filled = false;
		seqCount = 0;
		sumDust = maxDust >= 0 ? 0 : -1;
		
		if (seqMarks != null) {
			for (int i = 0; i < size; i++) {
				seqMarks[i] = 0;
			}
		}
	}

	public boolean isFilled() {
		return filled;
	}

	public int getSize() {
		return size;
	}

	public boolean isDust() {
		return sumDust > maxDust;
	}

	public int getDustValue() {
		return sumDust;
	}
	
	public int getMaxDust() {
		return maxDust;
	}
}