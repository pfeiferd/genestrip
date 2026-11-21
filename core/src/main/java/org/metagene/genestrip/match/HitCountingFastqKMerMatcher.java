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
package org.metagene.genestrip.match;

import org.metagene.genestrip.ExecutionContext;
import org.metagene.genestrip.store.KMerStore;
import org.metagene.genestrip.tax.SmallTaxTree;
import org.metagene.genestrip.tax.SmallTaxTree.SmallTaxIdNode;
import org.metagene.genestrip.util.LargeBitVector;

/**
 * A matcher that keeps a two-bit counter per stored k-mer instead of a single bit, so that the k-mers
 * hit exactly once and exactly twice can be counted per tax id ({@code kmersHitOnce},
 * {@code kmersHitTwice}). {@link SharedFractionEstimator} uses them for Chao's estimate of the shared
 * fraction, which, unlike the estimate from a single bit, is not thrown off by k-mers hit very often.
 * <p>
 * The counters live where the bits would: in the database's entries, or, with
 * {@link #isOwnUniqueKMerBits()}, in a vector of this matcher's own with two bits per stored k-mer. A
 * database whose values leave no room for the second bit is counted in such a vector too, with a
 * warning.
 */
public class HitCountingFastqKMerMatcher extends FastqKMerMatcher {
	/**
	 * Creates a matcher; the parameters are those of
	 * {@link FastqKMerMatcher#FastqKMerMatcher(KMerStore, int, int, ExecutionContext, boolean, SmallTaxTree, int, double, double, boolean, int, String)}.
	 *
	 * @param kmerStore              the k-mer store to match against
	 * @param initialReadSize        the initial read buffer size in bytes
	 * @param maxQueueSize           the maximum size of the read queue
	 * @param bundle                 the execution context
	 * @param withProbs              whether quality/probability information is processed
	 * @param taxTree                the taxonomy tree used for voting and LCA resolution
	 * @param maxPaths               the maximum number of candidate taxonomic paths per read
	 * @param maxReadTaxErrorCount   the maximum allowed tax error count per read
	 * @param maxReadClassErrorCount the maximum allowed classification error count per read
	 * @param writeAll               whether to write all reads, not just classified ones
	 * @param threshold              the minimum vote count threshold for classification
	 * @param dbMD5                  the MD5 checksum identifying the database
	 */
	public HitCountingFastqKMerMatcher(KMerStore<SmallTaxIdNode> kmerStore, int initialReadSize, int maxQueueSize,
			ExecutionContext bundle, boolean withProbs, SmallTaxTree taxTree, int maxPaths,
			double maxReadTaxErrorCount, double maxReadClassErrorCount, boolean writeAll, int threshold,
			String dbMD5) {
		super(kmerStore, initialReadSize, maxQueueSize, bundle, withProbs, taxTree, maxPaths, maxReadTaxErrorCount,
				maxReadClassErrorCount, writeAll, threshold, dbMD5);
	}

	@Override
	protected void initUniqueCounter() {
		if (!isOwnUniqueKMerBits()) {
			// Claims the database's marks, or fails as a plain matcher would.
			super.initUniqueCounter();
			if (kmerStore.setCountHits(true)) {
				return;
			}
			kmerStore.setMarkVisited(false);
			if (getLogger().isWarnEnabled()) {
				getLogger().warn("The database holds " + kmerStore.getNValues()
						+ " values, too many to count hits in its entries. Counting them in a vector of this"
						+ " matcher's own instead, at two bits per stored k-mer.");
			}
		}
		// Two bits per stored k-mer, counter i at bits 2i and 2i + 1 (see LargeBitVector#incrementCounter).
		if (uniqueKmerBits == null) {
			uniqueKmerBits = new LargeBitVector(2 * kmerStore.getEntries());
		} else {
			uniqueKmerBits.clear();
		}
	}

	@Override
	protected void markOwn(long pos) {
		uniqueKmerBits.incrementCounter(pos);
	}

	@Override
	protected void computeUniqueKmerCounts() {
		long[] hit = new long[statsIndex.length];
		long[] once = new long[statsIndex.length];
		long[] twice = new long[statsIndex.length];
		if (uniqueKmerBits != null) {
			final LargeBitVector counters = uniqueKmerBits;
			kmerStore.visit((store, kmer, index, pos) -> {
				final int counter = counters.getCounter(pos);
				if (counter != 0) {
					hit[index]++;
					if (counter == 1) {
						once[index]++;
					} else if (counter == 2) {
						twice[index]++;
					}
				}
			});
		} else {
			kmerStore.countHitsPerValueIndex(hit, once, twice);
		}
		for (int vi = 0; vi < statsIndex.length; vi++) {
			if (statsIndex[vi] != null) {
				statsIndex[vi].uniqueKmers = hit[vi];
				statsIndex[vi].kmersHitOnce = once[vi];
				statsIndex[vi].kmersHitTwice = twice[vi];
			}
		}
	}
}
