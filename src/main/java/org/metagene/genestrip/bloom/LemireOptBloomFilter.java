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
