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

public class CGATRingBufferTest extends TestCase {
	@Test
	public void testDust() {
		CGATRingBuffer buffer = new CGATRingBuffer(31, 20);
		int count = buffer.getSize();
		int[] fib = new int[count];
		if (count > 0) {
			fib[0] = 1;
		}
		if (count > 1) {
			fib[1] = 1;
		}
		for (int i = 2; i < fib.length; i++) {
			fib[i] = fib[i - 1] + fib[i - 2];
		}

		Random random = new Random(10);
		int dustCount = 0;
		int totalCount = 0;

		for (int j = 0; j < 1000; j++) {
			buffer.put((byte) 'N');
			assertEquals(0, buffer.getDustValue());
			for (int k = 0; k < 100000; k++) {
				assertEquals(k >= buffer.getSize(), buffer.isFilled());
				
				byte c = CGAT.DECODE_TABLE[random.nextInt(4)];
				buffer.put(c);

				if (buffer.isFilled()) {
					totalCount++;
					
					int n = 0;
					int sumDust = 0;
					int last = buffer.get(0);
					for (int i = 1; i < buffer.getSize(); i++) {
						byte d = buffer.get(i);
						if (d == last) {
							n++;
						} else if (n > 0) {
							sumDust += fib[n];
							n = 0;
						}
						last = d;
					}
					if (n > 0) {
						sumDust += fib[n];
					}

					assertEquals(buffer.getMaxDust() >= 0 ? sumDust : -1, buffer.getDustValue());
					/*
					if (sumDust > buffer.getMaxDust()) {
						System.out.println("Dust value: " + sumDust + " " + buffer.toString());
					}
					 */
					assertEquals(sumDust > buffer.getMaxDust() && buffer.getMaxDust() >= 0, buffer.isDust());
					if (buffer.isDust()) {
						dustCount++;
//						System.out.println(j + ": " + buffer + " " + sumDust);
					}

					// This just ensures that the dust measure is (really) symmetrical:
					int sumDust1 = sumDust;
					n = 0;
					sumDust = 0;
					last = buffer.get(buffer.getSize() - 1);
					for (int i = buffer.getSize() - 2; i >= 0; i--) {
						byte d = buffer.get(i);
						if (d == last) {
							n++;
						} else if (n > 0) {
							sumDust += fib[n];
							n = 0;
						}
						last = d;
					}
					if (n > 0) {
						sumDust += fib[n];
					}
					assertEquals(sumDust1, sumDust);
				}
			}
		}
		System.out.println("Dust ratio: " + ((double) dustCount) / totalCount);
	}

	public void testSomeKmerStrings() {
		CGATRingBuffer buffer = new CGATRingBuffer(31, 1000);
		fill(buffer, "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC");
		assertEquals(fib(31), buffer.getDustValue());
		fill(buffer, "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT");
		assertEquals(fib(31), buffer.getDustValue());
		fill(buffer, "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTA");
		assertEquals(fib(30), buffer.getDustValue());
		fill(buffer, "ATTTTTTTTTTTTTTTTTTTTTTTTTTTTTA");
		assertEquals(fib(29), buffer.getDustValue());
		fill(buffer, "TTTTTTTTTTTTTTTTTCCCCCCCCCCCCCC");
		assertEquals(fib(17) + fib(14), buffer.getDustValue());
		fill(buffer, "AAAATTTTTTTTTTTTTCCCCCCCCCCCCCC");
		assertEquals(fib(4) + fib(13) + fib(14), buffer.getDustValue());
		fill(buffer, "ACACACACACACACACACACACACACACACA");
		assertEquals(0, buffer.getDustValue());
	}

	public int fib(int n) {
		int f1 = 0, f2 = 1;
		int h;
		for (int i = 0; i < n; i++) {
			h = f2;
			f2 = f1 + f2;
			f1 = h;
		}
		return f1;
	}

	protected void fill(CGATRingBuffer buffer, String kmer) {
		for (int i = 0; i < kmer.length(); i++) {
			buffer.put((byte) kmer.charAt(i));
		}
	}

	@Test
	public void testGetKMer() {
		CGATRingBuffer buffer = new CGATRingBuffer(31, 25);
		CGATLongBuffer longBuffer = new CGATLongBuffer(31, 25);

		Random random = new Random(10);

		for (int i = 0; i < 1000; i++) {
			buffer.put((byte) 'N');
			longBuffer.put((byte) 'N');
			for (int j = 0; j < 1000; j++) {
				if (j < buffer.getSize()) {
					assertFalse(buffer.isFilled());
					assertEquals(-1, buffer.getKMer());

					assertFalse(longBuffer.isFilled());
					assertEquals(-1, longBuffer.getKMer());
				}
				byte c = CGAT.DECODE_TABLE[random.nextInt(4)];
				buffer.put(c);
				longBuffer.put(c);

				assertEquals(longBuffer.getKMer(), buffer.getKMer());
			}
		}
	}
}
