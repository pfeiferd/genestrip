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
		
		CGATRingBuffer buffer = new CGATRingBuffer(10, -1);
		
		Random random = new Random(10);
		
		for (int j = 0; j < 10000; j++) {
			byte c = CGAT.DECODE_TABLE[random.nextInt(4)];
			buffer.put(c);
			
			if (buffer.isFilled()) {
				int n = 0;
				int sumDust = 0;
				int last = buffer.get(0);
				if (j == 9998) {
					System.out.println(buffer);
				}
				for (int i = 1; i < buffer.getSize(); i++) {
					byte d = buffer.get(i);
					if (d == last) {
						n++;
					}
					else if (n > 0) {
						sumDust += fib[n];
						n = 0;
					}
					last = d;
				}
				if (n > 0) {
					sumDust += fib[n];
				}
				System.out.println(buffer + " " + sumDust);
//				assertEquals(sumDust, buffer.getDustValue());
			}
		}

	}
}
