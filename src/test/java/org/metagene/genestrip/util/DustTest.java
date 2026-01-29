package org.metagene.genestrip.util;

import junit.framework.TestCase;
import org.junit.Ignore;

import java.util.Random;

@Ignore
public class DustTest extends TestCase {
    private int k = 31;
    private final int[] fib;

    public DustTest() {
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

    public void testNaiveDust() {
        CGATRingBuffer buffer = new CGATRingBuffer(k, 500);
        CGATRingBuffer reverseBuffer = new CGATRingBuffer(k, 500);
        Random random = new Random(10);

        int[] pickThresholds = { 0, 50, 100, 150 };
        String[] picks = new String[pickThresholds.length];
        int[] picksDust = new int[pickThresholds.length];
        int[] countsBelowThresholds = new int[pickThresholds.length];
        for (int i = 0; i < picksDust.length; i++) {
            picksDust[i] = Integer.MAX_VALUE;
        }
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
                    for (int i = 0; i < pickThresholds.length; i++) {
                        if (res >= pickThresholds[i] && res < picksDust[i]) {
                            picks[i] = buffer.toString();
                            picksDust[i] = res;
                        }
                        if (res >= pickThresholds[i]) {
                            countsBelowThresholds[i]++;
                        }
                    }
                    counts++;

//                    System.out.println(buffer);
//                    System.out.println(res);
                }
            }
        }
        System.out.println("Total runs: " + counts);
        for (int i = 0; i < pickThresholds.length; i++) {
            System.out.println(pickThresholds[i] + " -> " + picksDust[i] + ":" + picks[i] );
            System.out.println(pickThresholds[i] + " counts below: " + countsBelowThresholds[i] + " ratio: " + ((double) countsBelowThresholds[i]) / counts);
        }
    }

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
                d += fib[srl1 / 2];
                srl1 = 0;
            }
            if (c == l3) {
                srl2++;
            }
            else {
                d += fib[srl2 / 3];
                srl2 = 0;
            }
            l3 = l2; l2 = l1; l1 = c;
        }
        d += fib[srl0] + fib[srl1 / 2] + fib[srl2 / 3];
        return d;
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
}
