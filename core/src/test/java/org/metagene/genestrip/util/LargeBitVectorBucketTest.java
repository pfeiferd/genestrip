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

/**
 * Checks that no bucket grows beyond 1 MB, so that G1 never allocates one as a humongous object.
 */
public class LargeBitVectorBucketTest {
	@Test
	public void testBucketShiftIsCapped() {
		// The 17.7 GB fill Bloom filter of hd-pseudomonas: about 2.2G words. Before the cap this
		// gave 2^22-word (32 MB) buckets.
		assertEquals(17, LargeBitVector.computeBucketShift(2_200_000_000L, LargeBitVector.DEFAULT_MIN_BUCKET_BITS));
		assertEquals(17, LargeBitVector.computeBucketShift(Long.MAX_VALUE >>> 6, LargeBitVector.DEFAULT_MIN_BUCKET_BITS));
	}

	@Test
	public void testSmallVectorsUnchanged() {
		// Below the cap the lock floor still decides: floor(log2(words)) - minBucketBits.
		assertEquals(10, LargeBitVector.computeBucketShift(1L << 19, 9));
		assertEquals(0, LargeBitVector.computeBucketShift(100, 9));
	}

	@Test
	public void testIndexingAcrossBuckets() {
		// 2^18 words: two full 2^17-word buckets at a floor of 1 bucket.
		LargeBitVector v = new LargeBitVector(64L << 18, 1);
		long last = (64L << 18) - 1;
		long boundary = 64L << 17;
		v.set(0);
		v.set(boundary - 1);
		v.set(boundary);
		v.set(last);
		assertEquals(true, v.get(0));
		assertEquals(true, v.get(boundary - 1));
		assertEquals(true, v.get(boundary));
		assertEquals(true, v.get(last));
		assertEquals(false, v.get(boundary + 1));
	}
}
