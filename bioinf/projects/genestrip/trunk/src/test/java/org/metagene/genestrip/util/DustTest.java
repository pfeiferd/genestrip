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
			if (j == 228) {
				System.out.println("stop");
			}
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
			}
		}
	}
}
