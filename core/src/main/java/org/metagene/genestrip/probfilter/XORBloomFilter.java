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
 * The {@code XOR} variant of {@link BloomFilter}: it hashes a k-mer by a bare exclusive or with the
 * hash function's factor, which is fast and apparently "good enough" for hashing k-mers, and pays for
 * it by reducing with a modulo instead of the multiply-shift of {@link BloomFilter#reduce(long)}.
 * <p>
 * This is the family where the trade pays most clearly: it evaluates the hash once per hash function -
 * seven times per lookup at a 1% target - so the cheaper hash outweighs the division the modulo costs.
 * Measured over 1e7 k-mers at the same false-positive rate, a lookup takes about 33 ns against
 * {@link BloomFilter}'s 39 ns, i.e. some 15% less. The one-hash families hash only once, so there
 * the margin is small and depends on the backing - see {@link XORBlockedBloomFilter}.
 */
public class XORBloomFilter extends BloomFilter {
    private static final long serialVersionUID = 1L;

    /**
     * Creates a filter targeting the given false-positive probability, sized for
     * {@code expectedInsertions} k-mers.
     *
     * @param fpp the target false-positive probability
     * @param expectedInsertions the expected number of k-mers to be inserted
     */
    public XORBloomFilter(double fpp, long expectedInsertions) {
        super(fpp, expectedInsertions);
    }

    /**
     * XORs the k-mer with the hash function's factor - no mixing whatsoever, which is what makes it
     * fast.
     * <p>
     * <strong>This relies on {@link BloomFilter#reduce(long)} being a modulo</strong>, i.e.
     * on a reduction that consumes every bit of the hash. Since a k-mer keeps its entropy in the low
     * bits and this hash leaves it there, any reduction driven by the high bits - such as a
     * multiply-shift - degrades this filter's false-positive rate by orders of magnitude, and through
     * the store's filter-based deduplication that costs k-mers. See {@code reduce} for the measurements.
     *
     * @param x the k-mer, encoded as a {@code long}, to hash
     * @param i the index of the hash function to apply
     * @return the {@code i}-th hash of the given k-mer
     */
    @Override
    protected final long hash(long x, final int i) {
        return hashFactors[i] ^ x;
    }

    /**
     * Reduces by a modulo rather than by the multiply-shift of
     * {@link BloomFilter#reduce(long)}, because {@link #hash(long, int)} does not mix.
     * <p>
     * <strong>The modulo must stay.</strong> Every multiply-shift is driven by the <em>high</em> bits of
     * its input, whereas a k-mer keeps its entropy in the <em>low</em> bits and this hash leaves it
     * there. This filter is the most sensitive of the package to that, because unlike
     * {@link BlockedBloomFilter} and {@link SingleWordBloomFilter} it derives nothing else from
     * the hash that would fold the low bits back in. Measured over 1e6 k-mers at a 1% target, letting it
     * inherit the multiply-shift yields 16.7% at {@code k=31} and 100% - every lookup a hit - at
     * {@code k=16}, against 1.02% and 1.00% with the modulo. That is not merely inaccurate:
     * {@code AbstractKMerStore} uses this filter to deduplicate while filling, so its false positives
     * make the store <em>drop k-mers</em>, which showed up as k-mers missing from a stored database. The
     * modulo consumes all of the hash's bits and keeps that intact.
     * <p>
     * Prefixing the multiply-shift with a mixing multiplication does repair the rate, but then costs two
     * multiplications per hash function and the whole point of this variant - a hash that does not mix -
     * is gone; that is what {@link MurmurBloomFilter} is for. Should this ever be revisited, measure
     * on keys with only {@code 2k} significant bits - random 64-bit keys hide the problem entirely.
     *
     * @param v the hash value to reduce
     * @return the bit index in {@code [0, bits)} for the given hash value
     */
    @Override
    protected final long reduce(final long v) {
        return Math.abs(v % bits);
    }
}
