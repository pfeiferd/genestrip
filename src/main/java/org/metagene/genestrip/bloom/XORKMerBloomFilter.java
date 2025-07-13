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
