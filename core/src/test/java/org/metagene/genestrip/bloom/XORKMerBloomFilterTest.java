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
package org.metagene.genestrip.bloom;

import static org.junit.Assert.assertTrue;

import org.junit.Test;

public class XORKMerBloomFilterTest extends KMerBloomFilterTest {
	@Override
	protected KMerProbFilter createFilter(long size, double fpp) {
		KMerProbFilter res = new XORKMerBloomFilter(fpp);
		res.clear();
		res.ensureExpectedSize(size, isTestLarge());
		return res;
	}

	@Override
	protected boolean isTestLarge() {
		return false;
	}

	// Regression: when the XOR hash returns exactly Long.MIN_VALUE, reduce() used to negate it first
	// (Math.abs(Long.MIN_VALUE) stays negative), producing a negative bit index and an out-of-bounds
	// access. This happens in practice once a large filter is hammered with enough hashes. The filter
	// is seeded deterministically, so craft an input whose hash 0 is exactly Long.MIN_VALUE.
	@Test
	public void testHashOfLongMinValueDoesNotOverflow() {
		XORKMerBloomFilter filter = new XORKMerBloomFilter(0.0001);
		filter.ensureExpectedSize(1000, false);
		// hash(x, 0) == hashFactors[0] ^ x, so this x maps to Long.MIN_VALUE.
		long x = filter.hashFactors[0] ^ Long.MIN_VALUE;
		filter.putLong(x); // must not throw ArrayIndexOutOfBoundsException
		assertTrue(filter.containsLong(x));
	}
}
