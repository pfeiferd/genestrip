package org.metagene.genestrip.finertree.bloom;

import org.metagene.genestrip.bloom.BlockedKMerBloomFilter;

public class KMerIndexBlockedBloomFilter extends BlockedKMerBloomFilter  implements KMerIndexProbFilter{
    public KMerIndexBlockedBloomFilter(long expectedInsertions) {
        super(expectedInsertions);
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
