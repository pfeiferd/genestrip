package org.metagene.genestrip.finertree.bloom;


public class KMerIndexFilterHelper {
    public static long combine(final long data, final int index) {
        return data ^ ((long) index) ^ (((long) index) << 32);
    }
}
