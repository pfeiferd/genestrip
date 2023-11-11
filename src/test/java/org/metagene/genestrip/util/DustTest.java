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

import static org.junit.Assert.assertEquals;

import java.util.Random;

import org.junit.Test;

public class DustTest {
	@Test
	public void testDustInCGATRingBuffer() {
		int count = 10;
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
		CGATRingBuffer buffer = new CGATRingBuffer(31, 25);

		Random random = new Random(10);

		for (int j = 0; j < 100000; j++) {
			byte c = CGAT.DECODE_TABLE[random.nextInt(4)];
			buffer.put(c);

			if (buffer.isFilled()) {
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
				assertEquals(sumDust > buffer.getMaxDust() && buffer.getMaxDust() >= 0, buffer.isDust());
				if (buffer.isDust()) {
					System.out.println(j + ": " + buffer + " " + sumDust);
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
}
