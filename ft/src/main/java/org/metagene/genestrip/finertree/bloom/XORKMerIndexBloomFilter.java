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

import org.metagene.genestrip.bloom.XORKMerBloomFilter;

/**
 * A XOR *k*-mer bloom filter that stores (k-mer, index) pairs rather than plain *k*-mers, by folding
 * the integer index into the *k*-mer value before hashing. Used to detect duplicate (k-mer, tax-id
 * index) combinations. The plain {@code putLong}/{@code containsLong} operations are therefore
 * unsupported.
 */
public class XORKMerIndexBloomFilter extends XORKMerBloomFilter {
    /**
     * Creates a filter sized for the given target false-positive probability.
     *
     * @param fpp the desired false-positive probability
     */
    public XORKMerIndexBloomFilter(double fpp) {
        super(fpp);
    }

    /**
     * @throws UnsupportedOperationException always; use {@link #putLongInt(long, int)} instead
     */
    @Override
    public void putLong(long data) {
        throw new UnsupportedOperationException();
    }

    /**
     * @throws UnsupportedOperationException always; use {@link #containsLongInt(long, int)} instead
     */
    @Override
    public boolean containsLong(long data) {
        throw new UnsupportedOperationException();
    }

    /**
     * Adds the given (k-mer, index) pair to the filter.
     *
     * @param data  the k-mer value
     * @param index the index folded into the k-mer value
     */
    public void putLongInt(long data, final int index) {
        data = data ^ ((long) index) ^ (((long) index) << 32);
        super.putLong(data);
    }

    /**
     * Tests whether the given (k-mer, index) pair was probably added to the filter.
     *
     * @param data  the k-mer value
     * @param index the index folded into the k-mer value
     * @return whether the given (k-mer, index) pair is (probably) contained in the filter
     */
    public boolean containsLongInt(long data, int index) {
        data = data ^ ((long) index) ^ (((long) index) << 32);
        return super.containsLong(data);
    }
}
