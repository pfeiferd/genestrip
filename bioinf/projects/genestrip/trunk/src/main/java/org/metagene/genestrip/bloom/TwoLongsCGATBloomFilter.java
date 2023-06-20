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
import java.util.ArrayList;
import java.util.List;

import org.metagene.genestrip.util.CGATRingBuffer;

public abstract class TwoLongsCGATBloomFilter extends AbstractCGATBloomFilter implements Serializable {
	private static final long serialVersionUID = 1L;

	private final int reverseEvenOffset;
	private final int reverseOddOffset;

//  Experimental stuff: minimize local dependencies of close k-mer bases with regards to generating hash code
//  by using this array that could help to access the bases in a non-sequential way.
//	protected final int[] longDistancePerm1;
//	protected final int[] longDistancePerm2;

	public TwoLongsCGATBloomFilter(int k, long expectedInsertions, double fpp) {
		super(k, expectedInsertions, fpp);

		reverseEvenOffset = k - 1;
		reverseOddOffset = k - 2;

		List<Integer> indices1 = new ArrayList<Integer>();
		List<Integer> indices2 = new ArrayList<Integer>();
		for (int i = 0; i < k; i++) {
			if (i % 2 == 0) {
				indices1.add(i);
			} else {
				indices2.add(i);
			}
		}

//		longDistancePerm1 = computBestPerm(indices1);
//		longDistancePerm2 = computBestPerm(indices2);
	}

//	private int[] computBestPerm(List<Integer> indices) {
//		double minValue = Double.MAX_VALUE;
//		int[] bestPerm = new int[indices.size()];
//		
//		int[] res = new int[indices.size()];
//		for (int index : indices) {
//			double v = computeMaxDistance(index, new ArrayList<Integer>(indices), res);
//			if (v < minValue) {
//				minValue = v;
//				System.arraycopy(res, 0, bestPerm, 0, bestPerm.length);
//			}
//		}
//		return res;
//	}
//
//	private double computeMaxDistance(int startIndex, List<Integer> indices, int[] res) {
//		res[0] = startIndex;
//		indices.remove((Object) startIndex);
//
//		double minValue = Double.MAX_VALUE;
//		int optimalIndex = 0;
//
//		for (int i = 1; i < res.length; i++) {
//			for (int index : indices) {
//				res[i] = index;
//				double d = computeDistance(res, i + 1);
//				if (d < minValue) {
//					minValue = d;
//					optimalIndex = index;
//				}
//			}
//			res[i] = optimalIndex;
//			indices.remove((Object) optimalIndex);
//		}
//
//		return minValue;
//	}
//
//	private double computeDistance(int[] solution, int size) {
//		double res = 1;
//		for (int i = 1; i < size; i++) {
//			res *= 1d / (solution[i - 1] - solution[i]);
//		}
//		return Math.abs(res);
//	}

	public void put(CGATRingBuffer buffer, int[] badPos) {
		long hash1 = hashStraight(buffer, true, badPos);
		if (hash1 == 0) {
			return;
		}
		long hash2 = hashStraight(buffer, false, badPos);
		if (hash2 == 0) {
			return;
		}
		putViaHash(hash1, hash2);
	}

	public void put(byte[] seq, int start, int[] badPos) {
		long hash1 = hashStraight(seq, start, true, badPos);
		if (hash1 == 0) {
			return;
		}
		long hash2 = hashStraight(seq, start, false, badPos);
		if (hash2 == 0) {
			return;
		}
		putViaHash(hash1, hash2);
	}

	protected long hashStraight(byte[] seq, int start, boolean even, int[] badPos) {
		long hash = 0;
		int c;

		int max = start + k;
		for (int i = start + (even ? 0 : 1); i < max; i += 2) {
			hash = hash << 2;
			c = CGAT_JUMP_TABLE[seq[i]];
			if (c == -1) {
				if (badPos != null) {
					badPos[0] = i;
				}
				return 0;
			}
			hash += c;
		}

		return hash == 0 ? 1L : hash;
	}

	protected long hashReverse(byte[] seq, int start, boolean even, int[] badPos) {
		long hash = 0;
		int c;

		for (int i = start + (even ? reverseEvenOffset : reverseOddOffset); i >= start; i -= 2) {
			hash = hash << 2;
			c = CGAT_REVERSE_JUMP_TABLE[seq[i]];
			if (c == -1) {
				if (badPos != null) {
					badPos[0] = i;
				}
				return 0;
			}
			hash += c;
		}

		return hash == 0 ? 1L : hash;
	}

	protected long hashStraight(CGATRingBuffer buffer, boolean even, int[] badPos) {
		long hash = 0;
		int c;

		for (int i = even ? 0 : 1; i < k; i += 2) {
			hash = hash << 2;
			c = CGAT_JUMP_TABLE[buffer.get(i)];
			if (c == -1) {
				if (badPos != null) {
					badPos[0] = i;
				}
				return 0;
			}
			hash += c;
		}

		return hash == 0 ? 1L : hash;
	}

	protected long hashReverse(CGATRingBuffer buffer, boolean even, int[] badPos) {
		long hash = 0;
		int c;

		for (int i = even ? reverseEvenOffset : reverseOddOffset; i >= 0; i -= 2) {
			hash = hash << 2;
			c = CGAT_REVERSE_JUMP_TABLE[buffer.get(i)];
			if (c == -1) {
				if (badPos != null) {
					badPos[0] = i;
				}
				return 0;
			}
			hash += c;
		}

		return hash == 0 ? 1L : hash;
	}

	protected long hash(CGATRingBuffer buffer, boolean even, boolean reverse, int[] badPos) {
		assert (buffer.getSize() == k);

		long hash = 0;
		int c;

		if (reverse) {
			for (int i = even ? 0 : 1; i < k; i += 2) {
				hash = hash << 2;
				c = CGAT_REVERSE_JUMP_TABLE[buffer.get(k - i - 1)];
				if (c == -1) {
					if (badPos != null) {
						badPos[0] = k - i - 1;
					}
					return 0;
				}
				hash += c;
			}
		} else {
			for (int i = even ? 0 : 1; i < k; i += 2) {
				hash = hash << 2;
				c = CGAT_JUMP_TABLE[buffer.get(i)];
				if (c == -1) {
					if (badPos != null) {
						badPos[0] = i;
					}
					return 0;
				}
				hash += c;
			}
		}
		return hash == 0 ? 1L : hash;
	}

	protected void putViaHash(long hash1, long hash2) {
		int hash;
		int index;

		for (int i = 0; i < hashes; i++) {
			hash = combineLongHashes(hash1, hash2, i);
			index = (int) ((hash >>> 6) % bits.length);
			bits[index] = bits[index] | (1L << (hash & 0b111111));
		}
	}

	protected abstract int combineLongHashes(long hash1, long hash2, int i);

	@Override
	public boolean containsStraight(byte[] seq, int start, int[] badPos) {
		long hash1 = hashStraight(seq, start, true, badPos);
		if (hash1 == 0) {
			return false;
		}
		long hash2 = hashStraight(seq, start, false, badPos);
		if (hash2 == 0) {
			return false;
		}
		return containsHash(hash1, hash2);
	}

	@Override
	public boolean containsReverse(byte[] seq, int start, int[] badPos) {
		long hash1 = hashReverse(seq, start, true, badPos);
		if (hash1 == 0) {
			return false;
		}
		long hash2 = hashReverse(seq, start, false, badPos);
		if (hash2 == 0) {
			return false;
		}
		return containsHash(hash1, hash2);
	}

	protected boolean containsHash(long hash1, long hash2) {
		int hash;
		int index;

		for (int i = 0; i < hashes; i++) {
			hash = combineLongHashes(hash1, hash2, i);
			index = (int) ((hash >>> 6) % bits.length);
			if (((bits[index] >> (hash & 0b111111)) & 1L) == 0) {
				return false;
			}
		}
		return true;
	}

	public boolean containsStraight(CGATRingBuffer buffer, int[] badPos) {
		long hash1 = hashStraight(buffer, true, badPos);
		if (hash1 == 0) {
			return false;
		}
		long hash2 = hashStraight(buffer, false, badPos);
		if (hash2 == 0) {
			return false;
		}
		return containsHash(hash1, hash2);
	}

	public boolean containsReverse(CGATRingBuffer buffer, int[] badPos) {
		long hash1 = hashReverse(buffer, true, badPos);
		if (hash1 == 0) {
			return false;
		}
		long hash2 = hashReverse(buffer, false, badPos);
		if (hash2 == 0) {
			return false;
		}
		return containsHash(hash1, hash2);
	}
}