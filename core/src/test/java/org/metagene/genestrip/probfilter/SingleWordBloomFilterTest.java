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
package org.metagene.genestrip.probfilter;

/**
 * Runs the shared {@link BloomFilterTest} suite against {@link SingleWordBloomFilter}, i.e.
 * checks that it never produces a false negative and stays within its false-positive budget. The
 * sizing is picked so that the filter meets the suite's 1% budget: confining a key's bits to one word
 * costs accuracy, so it needs a few bits per key more than {@link BlockedBloomFilter} for the
 * same rate (see {@link SingleWordBloomFilterInternalsTest} for the specifics).
 */
public class SingleWordBloomFilterTest extends BloomFilterTest {
	@Override
	protected double createFpp() {
		return 0.01;
	}

	@Override
	protected ProbFilter createFilter(long size, double fpp) {
		// The filter picks its backing from the size relative to MAX_SMALL_CAPACITY, which no
		// test-sized filter reaches, so ask for the bucketed backing explicitly.
		return isTestLarge() ? SingleWordBloomFilter.newLargeBacked(size, 16)
				: new SingleWordBloomFilter(size, 16);
	}

	@Override
	protected boolean isTestLarge() {
		return false;
	}
}
