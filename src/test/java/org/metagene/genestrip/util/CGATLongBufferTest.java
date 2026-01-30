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

import java.io.PrintStream;
import java.util.Random;

import org.junit.Ignore;
import org.junit.Test;

import junit.framework.TestCase;

public class CGATLongBufferTest extends TestCase {
	private int k = 31;
	private final int[] fib;

	public CGATLongBufferTest() {
		fib = new int[k];
		if (k > 0) {
			fib[0] = 0;
		}
		if (k > 1) {
			fib[1] = 1;
		}
		if (k > 2) {
			fib[2] = 2;
		}
		for (int i = 3; i < fib.length; i++) {
			fib[i] = fib[i - 1] + fib[i - 2];
		}
	}

	@Test
	public void testHighDustKMers() {
		CGATRingBuffer buffer = new CGATRingBuffer(31, 500);

		int v = fib[29];
		fill(buffer, 'A', 'C');
		assertEquals(v, buffer.getDustValue());
		fill(buffer, 'G', 'C');
		assertEquals(v, buffer.getDustValue());

		v = fib[30] + fib[29] + fib[28];
		fill(buffer, 'C', 'C');
		assertEquals(v, buffer.getDustValue());
		fill(buffer, 'T', 'T');
		assertEquals(v, buffer.getDustValue());
	}

	protected void fill(CGATRingBuffer buffer, char a, char b) {
		for (int i = 0; i < buffer.getSize(); i+=2) {
			buffer.putForTest((byte) a);
			buffer.putForTest((byte) b);
		}
		if (buffer.getSize() % 2 == 1) {
			buffer.putForTest((byte) a);
		}
	}

	@Test
	public void testSomeShortKmerStrings() {
		// From system documentation:
		CGATLongBuffer buffer = new CGATLongBuffer(8, 1000);
		fill(buffer, "TTTCGCGA");
		int v = fib[2] + fib[1] + fib[2];
		assertEquals(v, buffer.getDustValue());

		buffer = new CGATLongBuffer(4, 1000);
		fill(buffer, "ACAT");
		assertEquals(fib[1], buffer.getDustValue());
		fill(buffer, "AAAT");
		assertEquals( fib[2] + fib[1], buffer.getDustValue());
		fill(buffer, "AAAA");
		assertEquals(fib[3] + fib[2] + fib[1], buffer.getDustValue());
		fill(buffer, "ACAA");
		assertEquals(fib[1] + fib[1] + fib[1], buffer.getDustValue());

		buffer = new CGATLongBuffer(5, 1000);
		fill(buffer, "ACATA");
		assertEquals(fib[1] + fib[1], buffer.getDustValue());
		fill(buffer, "AAAAA");
		assertEquals(fib[4] + fib[3] + fib[2], buffer.getDustValue());
		fill(buffer, "AATAA");
		assertEquals(2 * fib[1] + fib[1] + fib[2], buffer.getDustValue());
		fill(buffer, "TATAT");
		assertEquals(fib[3], buffer.getDustValue());
	}

	protected void fill(CGATLongBuffer buffer, String kmer) {
		buffer.reset();
		for (int i = 0; i < kmer.length(); i++) {
			buffer.put((byte) kmer.charAt(i));
		}
	}

	@Test
	public void testViaNaiveDust() {
		CGATRingBuffer buffer = new CGATRingBuffer(k, 500);
		CGATRingBuffer reverseBuffer = new CGATRingBuffer(k, 500);
		Random random = new Random(10);

		int counts = 0;
		for (int j = 0; j < 1000; j++) {
			buffer.putForTest((byte) 'N');
			assertEquals(0, buffer.getDustValue());
			for (int h = 0; h < 100000; h++) {
				assertEquals(h >= buffer.getSize(), buffer.isFilled());

				byte c = CGAT.DECODE_TABLE[random.nextInt(4)];
				buffer.putForTest(c);

				if (buffer.isFilled()) {
					for (int i = 0; i < k; i++) {
						reverseBuffer.putForTest(buffer.get(k - i - 1));
					}
					int res = naiveDust(buffer);
					assertEquals(res, buffer.getDustValue());
					// Check for symmetry
					assertEquals(res, reverseBuffer.getDustValue());
				}
			}
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

	/* For Paper:
    // int fib[] = ... is assumed to be defined elsewhere.
    public int d(byte[] s) {
        int d = 0, srl0 = 0, srl1 = 0, srl2 = 0;
        int l1 = -1, l2 = -1, l3 = -1;
        for (int i = 0; i < s.length; i++) {
            if (s[i] == l1) srl0++;
            else {
                d += fib[srl0];
                srl0 = 0; }
            if (s[i] == l2) srl1++;
            else {
                d += fib[srl1 / 2];
                srl1 = 0; }
            if (s[i] == l3) srl2++;
            else {
                d += fib[srl2 / 3];
                srl2 = 0; }
            l3 = l2; l2 = l1; l1 = s[i];
        }
        d += fib[srl0] + fib[srl1 / 2] + fib[srl2 / 3];
        return d;
    }
    */

	@Ignore
	@Test
	public void testGenerateCSVTable() {
		int[] pickThresholds = new int[16];
		for (int i = 0; i < pickThresholds.length; i++) {
			pickThresholds[i] = fib[i];
		}
		String[] picks = new String[pickThresholds.length];
		int[] picksDust = new int[pickThresholds.length];
		int[] countsAboveThresholds = new int[pickThresholds.length];
		for (int i = 0; i < picksDust.length; i++) {
			picksDust[i] = Integer.MAX_VALUE;
		}

		CGATRingBuffer buffer = new CGATRingBuffer(k, 500);
		Random random = new Random(10);
		int max = 100000000;

		for (int j = 0; j < max; j++) {
			for (int h = 0; h < k; h++) {
				byte c = CGAT.DECODE_TABLE[random.nextInt(4)];
				buffer.putForTest(c);
			}
			int d = naiveDust(buffer);
			for (int i = 0; i < pickThresholds.length; i++) {
				if (d >= pickThresholds[i]) {
					countsAboveThresholds[i]++;
					if (d < picksDust[i]) {
						picks[i] = buffer.toString();
						picksDust[i] = d;
					}
				}
			}
		}
		PrintStream out = System.out;

		out.println("n; t = fib(n); counts >= t; >= t ratio; d(example); example");
		for (int i = 0; i < pickThresholds.length; i++) {
			out.print(i);
			out.print(";");
			out.print(pickThresholds[i]);
			out.print(";");
			out.print(countsAboveThresholds[i]);
			out.print(";");
			out.print(((double) countsAboveThresholds[i]) / max);
			out.print(";");
			out.print(picksDust[i]);
			out.print(";");
			out.print(picks[i]);
			out.println(";");
		}
	}

	// The non-streaming version of the dust function d - only for testing purposes.
	// Obviously to test the more complex streaming version CGATRingBuffer.
	protected int naiveDust(CGATRingBuffer buffer) {
		if (!buffer.isFilled()) {
			return -1;
		}
		int d = 0, srl0 = 0, srl1 = 0, srl2 = 0;
		int l1 = -1, l2 = -1, l3 = -1;
		for (int i = 0; i < buffer.getSize(); i++) {
			byte c = buffer.get(i);
			if (c == l1) {
				srl0++;
			}
			else {
				d += fib[srl0];
				srl0 = 0;
			}
			if (c == l2) {
				srl1++;
			}
			else {
				d += fib[srl1];
				srl1 = 0;
			}
			if (c == l3) {
				srl2++;
			}
			else {
				d += fib[srl2];
				srl2 = 0;
			}
			l3 = l2; l2 = l1; l1 = c;
		}
		d += fib[srl0] + fib[srl1] + fib[srl2];
		return d;
	}
}
