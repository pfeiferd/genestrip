/*
 *
 * "Commons Clause" License Condition v1.0
 *
 * The Software is provided to you by the Licensor under the License,
 * as defined below, subject to the following condition.
 *
 * Without limiting other conditions in the License, the grant of rights under the License
 * will not include, and the License does not grant to you, the right to Sell the Software.
 *
 * Software: genestrip
 * License: Apache 2.0
 * Licensor: Daniel Pfeifer (daniel.pfeifer@progotec.de)
 */
package org.metagene.genestrip.util;

import org.junit.Test;

import java.util.ArrayList;
import java.util.List;

import static org.junit.Assert.*;

public class SimpleBlockingQueueTest {

    // --- isEmpty() -----------------------------------------------------------

    @Test
    public void testIsEmptyOnNewQueue() {
        SimpleBlockingQueue<String> q = new SimpleBlockingQueue<>(4);
        assertTrue("New queue must be empty", q.isEmpty());
    }

    @Test
    public void testIsEmptyAfterPut() throws InterruptedException {
        SimpleBlockingQueue<String> q = new SimpleBlockingQueue<>(4);
        q.put("x");
        assertFalse("Queue with one element must not be empty", q.isEmpty());
    }

    @Test
    public void testIsEmptyAfterPutAndTake() throws InterruptedException {
        SimpleBlockingQueue<String> q = new SimpleBlockingQueue<>(4);
        q.put("x");
        q.take();
        assertTrue("Queue must be empty after taking all elements", q.isEmpty());
    }

    // --- size() / remainingCapacity() ----------------------------------------

    @Test
    public void testSizeAndRemainingCapacity() throws InterruptedException {
        int capacity = 5;
        SimpleBlockingQueue<Integer> q = new SimpleBlockingQueue<>(capacity);
        assertEquals(0, q.size());
        assertEquals(capacity, q.remainingCapacity());

        q.put(1);
        assertEquals(1, q.size());
        assertEquals(capacity - 1, q.remainingCapacity());

        q.put(2);
        assertEquals(2, q.size());

        q.take();
        assertEquals(1, q.size());
        assertEquals(capacity - 1, q.remainingCapacity());
    }

    // --- FIFO ordering -------------------------------------------------------

    @Test
    public void testFifoOrder() throws InterruptedException {
        SimpleBlockingQueue<Integer> q = new SimpleBlockingQueue<>(10);
        for (int i = 0; i < 10; i++) q.put(i);
        for (int i = 0; i < 10; i++)
            assertEquals("FIFO order violated at index " + i, Integer.valueOf(i), q.take());
    }

    // --- Circular wrap-around ------------------------------------------------

    @Test
    public void testCircularWrapAround() throws InterruptedException {
        SimpleBlockingQueue<Integer> q = new SimpleBlockingQueue<>(4);
        for (int i = 0; i < 4; i++) q.put(i);
        for (int i = 0; i < 4; i++) assertEquals(Integer.valueOf(i), q.take());
        for (int i = 10; i < 14; i++) q.put(i);
        for (int i = 10; i < 14; i++) assertEquals(Integer.valueOf(i), q.take());
        assertTrue(q.isEmpty());
    }

    // --- Multi-threaded producer / consumer -----------------------------------

    @Test
    public void testMultiThreadedProducerConsumer() throws InterruptedException {
        final int N = 1000;
        SimpleBlockingQueue<Integer> q = new SimpleBlockingQueue<>(16);
        List<Integer> received = new ArrayList<>(N);

        Thread producer = new Thread(() -> {
            for (int i = 0; i < N; i++) {
                try { q.put(i); } catch (InterruptedException e) { Thread.currentThread().interrupt(); }
            }
        });
        Thread consumer = new Thread(() -> {
            for (int i = 0; i < N; i++) {
                try { received.add(q.take()); } catch (InterruptedException e) { Thread.currentThread().interrupt(); }
            }
        });

        consumer.start();
        producer.start();
        producer.join(5000);
        consumer.join(5000);

        assertEquals("All items must be received", N, received.size());
        for (int i = 0; i < N; i++)
            assertEquals("Item out of order at " + i, Integer.valueOf(i), received.get(i));
    }

    // --- Blocking behaviour --------------------------------------------------

    @Test
    public void testTakeBlocksUntilPut() throws InterruptedException {
        SimpleBlockingQueue<String> q = new SimpleBlockingQueue<>(2);
        String[] result = new String[1];

        Thread consumer = new Thread(() -> {
            try { result[0] = q.take(); } catch (InterruptedException e) { Thread.currentThread().interrupt(); }
        });
        consumer.start();
        Thread.sleep(100);
        q.put("hello");
        consumer.join(2000);
        assertEquals("hello", result[0]);
    }

    @Test
    public void testPutBlocksWhenFull() throws InterruptedException {
        SimpleBlockingQueue<Integer> q = new SimpleBlockingQueue<>(2);
        q.put(1);
        q.put(2);

        boolean[] unblocked = {false};
        Thread producer = new Thread(() -> {
            try { q.put(3); unblocked[0] = true; }
            catch (InterruptedException e) { Thread.currentThread().interrupt(); }
        });
        producer.start();
        Thread.sleep(100);
        assertFalse("Producer should still be blocked", unblocked[0]);
        q.take();
        producer.join(2000);
        assertTrue("Producer must have continued after take()", unblocked[0]);
    }

    // --- Unsupported operations ----------------------------------------------

    @Test(expected = UnsupportedOperationException.class)
    public void testAddThrows() { new SimpleBlockingQueue<>(2).add("x"); }

    @Test(expected = UnsupportedOperationException.class)
    public void testOfferThrows() { new SimpleBlockingQueue<>(2).offer("x"); }

    @Test(expected = UnsupportedOperationException.class)
    public void testRemoveThrows() { new SimpleBlockingQueue<>(2).remove("x"); }

    @Test(expected = UnsupportedOperationException.class)
    public void testClearThrows() { new SimpleBlockingQueue<>(2).clear(); }

    @Test(expected = UnsupportedOperationException.class)
    public void testIteratorThrows() { new SimpleBlockingQueue<>(2).iterator(); }
}
