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

import org.junit.Test;

public class LargeBitVectorCounterTest {
	@Test
	public void testIncrementStopsAtThree() {
		LargeBitVector counters = new LargeBitVector(2 * 1000);
		for (int i = 0; i < 1000; i++) {
			for (int h = 0; h < i % 6; h++) {
				assertEquals(Math.min(h, 3), counters.incrementCounter(i));
			}
		}
		for (int i = 0; i < 1000; i++) {
			assertEquals("counter " + i, Math.min(i % 6, 3), counters.getCounter(i));
		}
		counters.clear();
		for (int i = 0; i < 1000; i++) {
			assertEquals(0, counters.getCounter(i));
		}
	}

	// Neighbouring counters share a word, and eight threads increment the same counters at once: a lost update
	// would leave a counter below its number of increments.
	@Test
	public void testConcurrentIncrements() throws InterruptedException {
		final int n = 100000;
		final LargeBitVector counters = new LargeBitVector(2L * n, 1);
		final int threads = 8;
		Thread[] workers = new Thread[threads];
		for (int t = 0; t < threads; t++) {
			workers[t] = new Thread(() -> {
				for (int i = 0; i < n; i += 3) {
					counters.incrementCounter(i);
				}
			});
			workers[t].start();
		}
		for (Thread worker : workers) {
			worker.join();
		}
		for (int i = 0; i < n; i++) {
			assertEquals("counter " + i, i % 3 == 0 ? 3 : 0, counters.getCounter(i));
		}
	}

	// Two increments per counter from two threads at once must reach exactly two.
	@Test
	public void testConcurrentTwo() throws InterruptedException {
		final int n = 200000;
		final LargeBitVector counters = new LargeBitVector(2L * n, 1);
		Thread[] workers = new Thread[2];
		for (int t = 0; t < 2; t++) {
			workers[t] = new Thread(() -> {
				for (int i = 0; i < n; i++) {
					counters.incrementCounter(i);
				}
			});
			workers[t].start();
		}
		for (Thread worker : workers) {
			worker.join();
		}
		for (int i = 0; i < n; i++) {
			assertEquals("counter " + i, 2, counters.getCounter(i));
		}
	}
}
