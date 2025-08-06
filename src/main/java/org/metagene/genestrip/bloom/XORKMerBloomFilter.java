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
    public XORKMerBloomFilter(int k, double fpp) {
        super(k, fpp);
    }

    @Override
    protected final long hash(long data, int i) {
        return hashFactors[i] ^ data;
    }

    // Optimization via inlining
    public final boolean containsLong(final long data) {
        for (int i = 0; i < hashes; i++) {
            // Inlined
            if (!bitVector.get(hashFactors[i] ^ data)) {
                return false;
            }
        }
        return true;
    }
}
