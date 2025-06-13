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

import java.util.Random;

import org.junit.Test;

import junit.framework.TestCase;

import static org.metagene.genestrip.util.CGAT.longToKMerStraight;

public class NextKMerTest extends TestCase {
	@Test
	public void testNextKMer() {
		byte[] kmerBuffer = new byte[31];
		for (int k = 1; k < kmerBuffer.length; k++) {
			CGATLongBuffer buffer = new CGATLongBuffer(k, -1);

			Random random = new Random(10);
			long oldKMer = -1;

			for (int j = 0; j < 100000; j++) {
				byte c = CGAT.DECODE_TABLE[random.nextInt(4)];
				buffer.put(c);

				if (buffer.isFilled()) {
					long kmer = buffer.getKMer();
					if (oldKMer != -1) {
						assertEquals(kmer, CGAT.nextKMerStraight(oldKMer, c, k));
					}
					oldKMer = kmer;
					longToKMerStraight(kmer, kmerBuffer, 0, k);
					String kmerStr = new String(kmerBuffer, 0, k);
					assertEquals(kmerStr, buffer.toString());
				}
			}
		}
		
		for (int k = 1; k < kmerBuffer.length; k++) {
			CGATLongBuffer buffer = new CGATLongBuffer(k, -1);

			Random random = new Random(10);
			long oldKMer = -1;

			for (int j = 0; j < 100000; j++) {
				byte c = CGAT.DECODE_TABLE[random.nextInt(4)];
				buffer.put(c);

				if (buffer.isFilled()) {
					long kmer = buffer.getReverseKMer();
					if (oldKMer != -1) {
						assertEquals(kmer, CGAT.nextKMerReverse(oldKMer, c, k));
					}
					oldKMer = kmer;
					longToKMerStraight(kmer, kmerBuffer, 0, k);
					CGAT.reverse(kmerBuffer,0, k);
					String kmerStr = new String(kmerBuffer, 0, k);
					assertEquals(kmerStr, buffer.toString());
				}
			}
		}
	}
	
	private static long kMerToLongStraight(CGATRingBuffer buffer) {
		long res = 0;
		int c;
		int k = buffer.getSize();
		for (int i = 0; i < k; i++) {
			res = Long.rotateLeft(res, 2);
			c = CGAT.CGAT_JUMP_TABLE[buffer.get(i)];
			res += c;
		}

		return res;
	}

	private static long kMerToLongReverse(CGATRingBuffer buffer) {
		long res = 0;
		int c;
		for (int i = buffer.getSize() - 1; i >= 0; i--) {
			res = Long.rotateLeft(res, 2);
			c = CGAT.CGAT_REVERSE_JUMP_TABLE[buffer.get(i)];
			res += c;
		}

		return res;
	}	
}
