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

// Added for expermimental reasons.
// This Bloom filter version requires a more complex hash function (from Lemire as below)
// than in the XORKMerBloomFilter which in combination makes it slower than the XORKMerBloomFilter.
public class LemireOptBloomFilter extends AbstractKMerBloomFilter {
    public LemireOptBloomFilter(double fpp) {
        super(fpp);
    }

    @Override
    protected final long hash(long x, final int i) {
        // return hashFactors[i] ^ x; -- makes test fail.
        x += hashFactors[i];
        x = (x ^ (x >>> 33)) * 0xff51afd7ed558ccdL;
        x = (x ^ (x >>> 33)) * 0xc4ceb9fe1a85ec53L;
        x = x ^ (x >>> 33);
        return x;
   }

    @Override
    protected long reduce(final long v) {
        if (largeBV) {
            return (v < 0 ? -v : v) % bits;
        } else {
            // Using optimization instead of '%', see:
            // http://lemire.me/blog/2016/06/27/a-fast-alternative-to-the-modulo-reduction/
            // and Line 34 in https://github.com/FastFilter/fastfilter_java/blob/master/fastfilter/src/main/java/org/fastfilter/utils/Hash.java
            // In general, it would NOT work for long indexes because of overflow due to the multiplicaton.
            return (int) (((((int) v) & 0xffffffffL) * (bits & 0xffffffffL)) >>> 32);
        }
    }
}
