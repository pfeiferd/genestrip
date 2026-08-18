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

/**
 * Selects one k-mer in n, as a function of the k-mer and of nothing else.
 * <p>
 * Wherever a fraction of the k-mers is enough - a database that is meant to hold only a share of
 * them, a count that is only estimated - the share has to be chosen by the k-mer itself. A k-mer is
 * then either selected in every genome it occurs in or in none of them, and that is what the two
 * users of this class rest on:
 * <ul>
 * <li>The database keeps one k-mer in {@code kMerSampling}. A k-mer's tax id is the lowest common
 * ancestor of the taxa of all regions it was met in, which the update pass computes by reading every
 * region of every taxon. Selecting by position instead - every n-th base of a region, as this was
 * done before - would meet a stored k-mer again only when it happened to fall on a selected offset
 * in the other genome, which for one in n fails with probability 1-1/n per occurrence. The stored
 * k-mers would come out looking more specific than they are, and reads of a taxon whose occurrence
 * was missed would be attributed to another one.</li>
 * <li>{@code kmerindexsize} estimates how many (k-mer, leaf) pairs the k-mer index filter will hold
 * by counting those of a sample. Counting the distinct pairs of the sample and scaling by the rate
 * only works if a whole k-mer is in or out: the pairs are partitioned by their k-mer, and a rule
 * that let a k-mer in sometimes would count occurrences rather than distinct pairs. Measured against
 * real RefSeq sequence at a duplication factor of 9, selecting every n-th pair by position
 * overestimated by 94 per cent at one in 4 and by 398 per cent at one in 64, where this rule stayed
 * within half a per cent.</li>
 * </ul>
 * <p>
 * Composing the two is exact rather than approximate, and that is deliberate. Both use this same
 * rule, so the finer sample is a <em>subset</em> of the coarser one and not an independent share of
 * it: of a population already reduced to one in {@code s}, a threshold for one in {@code q} keeps a
 * fraction of {@code s/max(s, q)}, not {@code 1/q}. A caller that samples what is already sampled
 * has to scale by that - see {@code KMerIndexSizeGoal} - which needs no assumption about two hash
 * functions being independent of one another.
 * <p>
 * The decision is a multiplication and a comparison, the whole of it, because it is put for every
 * k-mer of every reference sequence. Deciding on {@code kmer % n == 0} instead would be cheaper
 * still and is sound for odd {@code n}, but fails badly for even {@code n}: a k-mer is two bits per
 * base, so {@code kmer % 16} is its last two bases, and selecting the k-mers ending in {@code AA}
 * took 20.6 per cent of real RefSeq sequence rather than 6.25 and overestimated a count by 204 per
 * cent. Running a MurmurHash3 finalizer is equally accurate but twice as expensive (4.9 against 2.4
 * ns per k-mer). The product's upper bits depend on every bit of the k-mer, which is all a yes/no
 * decision needs; measured over 60 million real canonical k-mers the share selected was within 0.1
 * per cent of the rate for one in 16, 17 and 64 alike.
 * <p>
 * The k-mer handed in must be the canonical one (see
 * {@link CGATLongBuffer#getStandardKMer()}), since that is what a database stores: selecting on the
 * forward encoding would treat a sequence and its reverse complement as two different k-mers.
 */
public class KMerSampling {
    /**
     * The threshold that selects every k-mer, i.e. no sampling at all. It is not a value the
     * comparison could produce, so it doubles as the marker for the unsampled case and keeps that case
     * down to a single comparison.
     */
    public static final long ALL = -1L;

    /**
     * The odd constant the k-mer is multiplied by. Any odd constant makes the multiplication a
     * bijection; this one - 2^64 divided by the golden ratio - spreads the upper bits well, which is
     * the half of the product the threshold looks at.
     */
    private static final long MULTIPLIER = 0x9E3779B97F4A7C15L;

    private KMerSampling() {
    }

    /**
     * Returns the threshold that selects one k-mer in {@code oneIn}.
     *
     * @param oneIn one in how many k-mers to select; 1 or less selects every one of them
     * @return the threshold to pass to {@link #isSampled(long, long)}
     */
    public static long thresholdForOneIn(int oneIn) {
        return oneIn <= 1 ? ALL : Long.divideUnsigned(-1L, oneIn);
    }

    /**
     * Returns whether the given k-mer belongs to the sample described by the threshold.
     * <p>
     * Callers hold the threshold in a field of their own and pass it in, so that this stays a static
     * call over two values on the hot path.
     *
     * @param kmer      the canonical k-mer, encoded as a {@code long}
     * @param threshold a threshold from {@link #thresholdForOneIn(int)}, or {@link #ALL}
     * @return whether the k-mer is part of the sample
     */
    public static boolean isSampled(long kmer, long threshold) {
        return threshold == ALL || Long.compareUnsigned(kmer * MULTIPLIER, threshold) < 0;
    }
}
