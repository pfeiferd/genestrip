package org.metagene.genestrip.finertree.bloom;

public interface KMerIndexProbFilter {
    /**
     * Adds the given (k-mer, index) pair to the filter and reports whether it was newly added, using
     * the same XOR-folding of the index into the k-mer as {@link #putLongInt(long, int)} but
     * computing each hash only once and setting the bits atomically. This combines a {@link
     * #containsLongInt(long, int)} check with the insertion in a single hashing pass and is safe for
     * concurrent use by multiple threads. The resulting filter state is identical to a plain {@code if
     * (!containsLongInt(data, index)) putLongInt(data, index)} sequence.
     *
     * @param data  the k-mer encoded as a long
     * @param index the integer index to associate with the k-mer
     * @return {@code true} if the (k-mer, index) pair was not already present, {@code false} otherwise
     */
    public boolean putLongInt(long data, final int index);

    /**
     * Tests whether the given (k-mer, index) pair was likely added to the filter, using the same
     * XOR-folding of the index into the k-mer as {@link #putLongInt(long, int)}. As with any Bloom
     * filter, false positives are possible but false negatives are not.
     *
     * @param data  the k-mer encoded as a long
     * @param index the integer index associated with the k-mer
     * @return {@code true} if the pair is possibly present, {@code false} if it is definitely absent
     */
    public boolean containsLongInt(long data, int index);

    public long getBitSize();
}
