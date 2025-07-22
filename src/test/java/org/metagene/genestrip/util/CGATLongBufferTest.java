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

public class CGATLongBufferTest extends TestCase {
	public void testHighDustKMers() {
		CGATLongBuffer buffer = new CGATLongBuffer(31, 500);

		int v = 832040;
		fill(buffer, 'A', 'C');
		assertEquals(v, buffer.getDustValue());
		fill(buffer, 'C', 'C');
		assertEquals(v, buffer.getDustValue());
		fill(buffer, 'G', 'C');
		assertEquals(v, buffer.getDustValue());
		fill(buffer, 'T', 'T');
	}

	protected void fill(CGATLongBuffer buffer, char a, char b) {
		for (int i = 0; i < buffer.getSize(); i+=2) {
			buffer.put((byte) a);
			buffer.put((byte) b);
		}
		if (buffer.getSize() % 2 == 1) {
			buffer.put((byte) a);
		}
	}



	@Test
	public void testDust() {
		CGATRingBuffer buffer = new CGATRingBuffer(31, 200);
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
		fib[0] = 0;

		Random random = new Random(10);
		int dustCount = 0;
		int totalCount = 0;

		for (int j = 0; j < 1000; j++) {
			buffer.putForTest((byte) 'N');
			assertEquals(0, buffer.getDustValue());
			for (int k = 0; k < 100000; k++) {
				assertEquals(k >= buffer.getSize(), buffer.isFilled());
				
				byte c = CGAT.DECODE_TABLE[random.nextInt(4)];
				buffer.putForTest(c);

				if (buffer.isFilled()) {
					totalCount++;
					
					int n = 0;
					int sumDust = 0;
					int beforeLast = -1;
					int last = -1;
					for (int i = 0; i < buffer.getSize(); i++) {
						byte d = buffer.get(i);
						if (d == beforeLast) {
							n++;
						} else {
							sumDust += fib[n];
							n = 0;
						}
						beforeLast = last;
						last = d;
					}
					sumDust += fib[n];
//					System.out.println(k + ": " + buffer);
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
					beforeLast = -1;
					last = buffer.get(buffer.getSize() - 1);
					for (int i = buffer.getSize() - 2; i >= 0; i--) {
						byte d = buffer.get(i);
						if (d == beforeLast) {
							n++;
						} else if (n > 0) {
							sumDust += fib[n];
							n = 0;
						}
						beforeLast = last;
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

	public void testSomeShortKmerStrings() {
		CGATLongBuffer buffer = new CGATLongBuffer(4, 1000);
		fill(buffer, "ACAT");
		assertEquals(fib(2), buffer.getDustValue());
		fill(buffer, "AAAT");
		assertEquals(fib(2), buffer.getDustValue());
		fill(buffer, "AAAA");
		assertEquals(fib(3), buffer.getDustValue());
		fill(buffer, "ACAA");
		assertEquals(fib(2), buffer.getDustValue());

		buffer = new CGATLongBuffer(5, 1000);
		fill(buffer, "ACATA");
		assertEquals(2* fib(2), buffer.getDustValue());
		fill(buffer, "AAAAA");
		assertEquals(fib(4), buffer.getDustValue());
		fill(buffer, "AATAA");
		assertEquals(fib(2), buffer.getDustValue());
		fill(buffer, "TATAT");
		assertEquals(fib(4), buffer.getDustValue());
	}

	public void testSomeKmerStrings() {
		CGATLongBuffer buffer = new CGATLongBuffer(31, 1000);
		fill(buffer, "ACACACACACACACACACACACACACACACA");
		assertEquals(fib(30), buffer.getDustValue());
		fill(buffer, "ATACACACACACACACACACACACACACACA");
		assertEquals(fib(1) + fib(28), buffer.getDustValue());
		fill(buffer, "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC");
		assertEquals(fib(30), buffer.getDustValue());
		fill(buffer, "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT");
		assertEquals(fib(30), buffer.getDustValue());
		fill(buffer, "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTA");
		assertEquals(fib(29), buffer.getDustValue());
		fill(buffer, "ATTTTTTTTTTTTTTTTTTTTTTTTTTTTTA");
		assertEquals(fib(28), buffer.getDustValue());
		fill(buffer, "TTTTTTTTTTTTTTTTTCCCCCCCCCCCCCC");
		assertEquals(fib(16) + fib(13), buffer.getDustValue());
		fill(buffer, "AAAATTTTTTTTTTTTTCCCCCCCCCCCCCC");
		assertEquals(fib(3) + fib(12) + fib(13), buffer.getDustValue());
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

	protected void fill(CGATLongBuffer buffer, String kmer) {
		buffer.reset();
		for (int i = 0; i < kmer.length(); i++) {
			buffer.put((byte) kmer.charAt(i));
		}
	}

	@Test
	public void testGetKMer() {
		CGATLongBuffer longBuffer = new CGATLongBuffer(31, 25);
		CGATRingBuffer ringBuffer = new CGATRingBuffer(31, 25);

		Random random = new Random(10);

		for (int i = 0; i < 1000; i++) {
			longBuffer.put((byte) 'N');
			ringBuffer.putForTest((byte) 'N');
			for (int j = 0; j < 1000; j++) {
				if (j < longBuffer.getSize()) {
					assertFalse(longBuffer.isFilled());
					assertEquals(-1, longBuffer.getKMer());
				}
				byte c = CGAT.DECODE_TABLE[random.nextInt(4)];
				longBuffer.put(c);
				ringBuffer.putForTest(c);

				assertEquals(longBuffer.getKMer(), ringBuffer.getKMerForTest());
			}
		}
	}
}
