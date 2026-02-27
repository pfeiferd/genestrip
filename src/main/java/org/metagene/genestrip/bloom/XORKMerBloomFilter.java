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

    public XORKMerBloomFilter(double fpp) {
        super(fpp);
    }

    @Override
    protected final long hash(long x, final int i) {
        return hashFactors[i] ^ x;
    }

    // The super simple hash function from above cannot be combined with Lemire's optimization.
    // It results in a very high fpp...
    // The more complex hash function from Lemire outweighs the cost of the modulo operator used here.
    // That's why we keep it simple and leave it as it is.
    protected final long reduce(final long v) {
        return (v < 0 ? -v : v) % bits;
    }
}
