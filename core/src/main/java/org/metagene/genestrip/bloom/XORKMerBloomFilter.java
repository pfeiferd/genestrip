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

/**
 * Fast and apparently "good enough" for hashing k-mers ...
 */
public class XORKMerBloomFilter extends AbstractKMerBloomFilter {
    private static final long serialVersionUID = 1L;

    /**
     * Creates a filter targeting the given false-positive probability, sized for
     * {@code expectedInsertions} k-mers.
     *
     * @param fpp the target false-positive probability
     * @param expectedInsertions the expected number of k-mers to be inserted
     */
    public XORKMerBloomFilter(double fpp, long expectedInsertions) {
        super(fpp, expectedInsertions);
    }

    /**
     * XORs the k-mer with the hash function's factor - no mixing whatsoever, which is what makes it
     * fast.
     * <p>
     * <strong>This relies on {@link AbstractKMerBloomFilter#reduce(long)} being a modulo</strong>, i.e.
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
}
