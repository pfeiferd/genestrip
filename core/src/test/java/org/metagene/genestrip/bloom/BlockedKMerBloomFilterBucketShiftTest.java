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

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.util.Random;

import org.junit.Test;

/**
 * Verifies the configurable large-backing bucket width of {@link BlockedKMerBloomFilter}: a small
 * bucket shift (so a key's block frequently straddles a bucket boundary) yields no false negatives, the
 * width survives serialization, and an out-of-range width is rejected.
 */
public class BlockedKMerBloomFilterBucketShiftTest {

	private static long[] randomKMers(int n, long seed) {
		Random random = new Random(seed);
		long[] kmers = new long[n];
		for (int i = 0; i < n; i++) {
			kmers[i] = random.nextLong();
		}
		return kmers;
	}

	/** Fills a large-backed filter of the given bucket shift with the keys and returns it. */
	private static BlockedKMerBloomFilter filledLargeFilter(int bucketShift, long[] kmers) {
		// Ask for the bucketed backing explicitly - these small test sizes would never reach it
		// through the size alone.
		BlockedKMerBloomFilter filter = BlockedKMerBloomFilter.newLargeBacked(kmers.length, 10, 42L, bucketShift);
		// Guards the seam itself: without this, a newLargeBacked() that quietly fell back to the small
		// backing would leave every test below passing while testing the wrong path.
		assertTrue("must use the bucketed backing", filter.isLargeBacked());
		for (long kmer : kmers) {
			filter.putLong(kmer);
		}
		return filter;
	}

	private static void assertNoFalseNegatives(BlockedKMerBloomFilter filter, long[] kmers) {
		for (long kmer : kmers) {
			assertTrue("inserted key must be found", filter.containsLong(kmer));
		}
	}

	private static byte[] serialize(BlockedKMerBloomFilter filter) throws IOException {
		ByteArrayOutputStream bytes = new ByteArrayOutputStream();
		try (ObjectOutputStream out = new ObjectOutputStream(bytes)) {
			out.writeObject(filter);
		}
		return bytes.toByteArray();
	}

	private static BlockedKMerBloomFilter deserialize(byte[] data) throws IOException, ClassNotFoundException {
		try (ObjectInputStream in = new ObjectInputStream(new ByteArrayInputStream(data))) {
			return (BlockedKMerBloomFilter) in.readObject();
		}
	}

	@Test
	public void testSmallBucketShiftHasNoFalseNegatives() {
		// A 2^10-word bucket over enough keys spans many buckets, so many keys' two-word blocks straddle a
		// bucket boundary - which must still be addressed correctly (no lost bits).
		long[] kmers = randomKMers(60_000, 7);
		assertNoFalseNegatives(filledLargeFilter(BlockedKMerBloomFilter.MIN_BUCKET_SHIFT, kmers), kmers);
	}

	@Test
	public void testCustomBucketShiftSurvivesSerialization() throws IOException, ClassNotFoundException {
		long[] kmers = randomKMers(60_000, 11);
		BlockedKMerBloomFilter loaded = deserialize(serialize(filledLargeFilter(12, kmers)));
		assertNoFalseNegatives(loaded, kmers);
	}

	@Test
	public void testBucketShiftOutOfRangeRejected() {
		try {
			new BlockedKMerBloomFilter(1000, 10, 42L, BlockedKMerBloomFilter.MAX_BUCKET_SHIFT + 1);
			fail("expected IllegalArgumentException for out-of-range bucketShift");
		} catch (IllegalArgumentException expected) {
			// expected
		}
	}

	@Test
	public void testMinBucketShiftUsesFewestBuckets() {
		// A sizing whose word count fits within one int-addressable block gets a single bucket: the
		// shift is the smallest power of two that still holds every word.
		int bitsPerKey = 10;
		for (long expected : new long[] { 1, 100, 100_000, 10_000_000 }) {
			int shift = BlockedKMerBloomFilter.minBucketShift(expected, bitsPerKey);
			long words = (Math.max(1, expected) * bitsPerKey + 63) / 64 + 16 + 1;
			assertTrue("one bucket must hold all words", (1L << shift) >= words);
			assertTrue("no smaller power of two would suffice",
					shift == BlockedKMerBloomFilter.MIN_BUCKET_SHIFT || (1L << (shift - 1)) < words);
		}
	}

	@Test
	public void testMinBucketShiftCapsAtMax() {
		// A sizing too large for a single int-addressable bucket falls back to the widest permitted
		// bucket, so the grid still uses as few buckets as the cap allows.
		long huge = (1L << 40) / 10; // words far exceed 2^MAX_BUCKET_SHIFT
		assertEquals(BlockedKMerBloomFilter.MAX_BUCKET_SHIFT,
				BlockedKMerBloomFilter.minBucketShift(huge, 10));
	}
}
