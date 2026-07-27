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
package org.metagene.genestrip.finertree.bloom;

import org.metagene.genestrip.bloom.KMerProbFilter;
import org.metagene.genestrip.bloom.XORKMerBloomFilter;

/**
 * A {@link XORKMerBloomFilter} that stores k-mers together with an associated integer index by
 * folding the index into the k-mer via XOR before delegating to the underlying filter. This lets a
 * single Bloom filter answer membership questions for (k-mer, index) pairs rather than plain
 * k-mers. The plain {@link #putLong(long)} and {@link #containsLong(long)} operations are therefore
 * disabled; callers must use {@link #putLongInt(long, int)} and {@link #containsLongInt(long, int)}.
 */
public class XORKMerIndexBloomFilter extends XORKMerBloomFilter implements KMerIndexProbFilter {
    /**
     * Creates an index Bloom filter with the given target false-positive probability, sized for
     * {@code expectedInsertions} (k-mer, index) pairs.
     *
     * @param fpp the desired false-positive probability of the filter
     * @param expectedInsertions the expected number of (k-mer, index) pairs to be inserted
     */
    public XORKMerIndexBloomFilter(double fpp, long expectedInsertions) {
        super(fpp, expectedInsertions);
    }

    /**
     * Unsupported: a plain k-mer cannot be queried because this filter always folds an index into
     * the k-mer. Use {@link #containsLongInt(long, int)} instead.
     *
     * @param data the k-mer encoded as a long
     * @return never returns normally
     * @throws UnsupportedOperationException always
     */
    @Override
    public boolean containsLong(long data) {
        throw new UnsupportedOperationException();
    }

    /**
     * Unsupported: a plain k-mer cannot be added because this filter always folds an index into the
     * k-mer. Use {@link #putLongInt(long, int)} instead.
     *
     * @param data the k-mer encoded as a long
     * @return never returns normally
     * @throws UnsupportedOperationException always
     */
    @Override
    public boolean putLong(long data) {
        throw new UnsupportedOperationException();
    }

    @Override
    public boolean putLongInt(long data, final int index) {
        data = data ^ ((long) index) ^ (((long) index) << 32);
        return super.putLong(data);
    }

    @Override
    public boolean containsLongInt(long data, int index) {
        data = data ^ ((long) index) ^ (((long) index) << 32);
        return super.containsLong(data);
    }
}
