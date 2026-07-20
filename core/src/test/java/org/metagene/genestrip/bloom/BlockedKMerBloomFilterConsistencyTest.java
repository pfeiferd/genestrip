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
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.util.Random;

import org.junit.Test;

/**
 * Verifies that {@link BlockedKMerBloomFilter#putLong(long)} — the combined single-pass insert
 * — behaves exactly like the previous {@code if (!containsLong(x)) putLong(x)} sequence (this filter is
 * filled single-threaded only), and that membership survives serialization.
 */
public class BlockedKMerBloomFilterConsistencyTest {

	private static BlockedKMerBloomFilter newFilter(long expected) {
		return newFilter(expected, false);
	}

	private static BlockedKMerBloomFilter newFilter(long expected, boolean large) {
		// 'large' forces the bucketed backing at small sizes by lowering the small/large threshold
		// the constructor consults (it is otherwise chosen from the size).
		long saved = BlockedKMerBloomFilter.MAX_SMALL_CAPACITY;
		if (large) {
			BlockedKMerBloomFilter.MAX_SMALL_CAPACITY = 1;
		}
		try {
			return new BlockedKMerBloomFilter(expected);
		} finally {
			BlockedKMerBloomFilter.MAX_SMALL_CAPACITY = saved;
		}
	}

	private static long[] randomKMers(int n, long seed) {
		Random random = new Random(seed);
		long[] kmers = new long[n];
		for (int i = 0; i < n; i++) {
			kmers[i] = random.nextLong();
		}
		return kmers;
	}

	/** Asserts the two filters answer every membership query identically over a large probe set. */
	private static void assertEquivalent(BlockedKMerBloomFilter expected, BlockedKMerBloomFilter actual,
			long[] insertedKeys, long probeSeed) {
		for (long kmer : insertedKeys) {
			assertEquals("membership on inserted key", expected.containsLong(kmer), actual.containsLong(kmer));
		}
		Random probe = new Random(probeSeed);
		for (int i = 0; i < insertedKeys.length; i++) {
			long randomKmer = probe.nextLong();
			assertEquals("membership on random key", expected.containsLong(randomKmer), actual.containsLong(randomKmer));
		}
	}

	@Test
	public void testMatchesClassicPathSmallBacking() {
		checkMatchesClassicPath(false);
	}

	@Test
	public void testMatchesClassicPathLargeBacking() {
		checkMatchesClassicPath(true);
	}

	private void checkMatchesClassicPath(boolean large) {
		int n = 20_000;
		long[] kmers = randomKMers(n, 314);

		BlockedKMerBloomFilter reference = newFilter(n, large);
		for (long kmer : kmers) {
			// Feed each k-mer twice to exercise the "already present" branch.
			if (!reference.containsLong(kmer)) {
				reference.putLong(kmer);
			}
			if (!reference.containsLong(kmer)) {
				reference.putLong(kmer);
			}
		}

		BlockedKMerBloomFilter candidate = newFilter(n, large);
		for (long kmer : kmers) {
			boolean firstAdded = candidate.putLong(kmer);
			boolean secondAdded = candidate.putLong(kmer);
			assertFalse("repeat insert must not be reported as new", secondAdded);
			if (!firstAdded) {
				assertTrue(candidate.containsLong(kmer));
			}
		}

		// The two backing words a key touches are always distinct, so single-threaded "newly added"
		// detection is exact and every membership answer must match.
		assertEquivalent(reference, candidate, kmers, 999);
	}

	@Test
	public void testSerializationPreservesMembership() throws IOException, ClassNotFoundException {
		int n = 5_000;
		long[] kmers = randomKMers(n, 4242);
		BlockedKMerBloomFilter filter = newFilter(n);
		for (long kmer : kmers) {
			filter.putLong(kmer);
		}

		ByteArrayOutputStream bytes = new ByteArrayOutputStream();
		try (ObjectOutputStream out = new ObjectOutputStream(bytes)) {
			out.writeObject(filter);
		}
		BlockedKMerBloomFilter loaded;
		try (ObjectInputStream in = new ObjectInputStream(new ByteArrayInputStream(bytes.toByteArray()))) {
			loaded = (BlockedKMerBloomFilter) in.readObject();
		}

		for (long kmer : kmers) {
			assertTrue("no false negative after reload", loaded.containsLong(kmer));
		}
	}
}
