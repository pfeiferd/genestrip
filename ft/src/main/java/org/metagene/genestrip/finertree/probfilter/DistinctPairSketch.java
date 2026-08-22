package org.metagene.genestrip.finertree.probfilter;

import net.agkn.hll.HLL;
import org.metagene.genestrip.util.MurmurHash3DropIn;

import java.util.Collection;
import java.util.concurrent.ConcurrentLinkedQueue;

/**
 * Counts how many distinct (<em>k</em>-mer, node) pairs a pass over the reference sequences produces,
 * by sketching them with HyperLogLog instead of storing them.
 * <p>
 * Two goals need this number and neither can afford to hold the pairs: {@code kmerindexsize} sizes
 * the filter of {@code kmerindexbloom} with it, and {@code dbqualcounts} sizes the filter it
 * deduplicates its counts through. Both otherwise fall back on a bound that assumes every {@code
 * k}-mer of a taxon to occur in every genome below it, which can exceed the truth by orders of
 * magnitude -- for a database whose genus holds hundreds of millions of {@code k}-mers over hundreds
 * of leaves, by enough to make the filter itself unallocatable.
 * <p>
 * The two goals read the sequences through different reader hierarchies and mean different things by
 * the second component of a pair, so they share this by holding one rather than by inheriting from a
 * common goal. What they do share is the whole of the counting: the sketch per reader thread, the
 * hashing, the merge and the scaling of a sample back up to the population.
 */
public class DistinctPairSketch {
    /** Registers of the sketches: 2^15 of 6 bits, for a relative error of about half a percent. */
    private static final int HLL_LOG2M = 15;
    private static final int HLL_REGISTER_WIDTH = 6;
    /**
     * The seed all sketches hash with. Merging two sketches is only meaningful when they hashed their
     * input alike, so this is a constant and not drawn per reader.
     */
    private static final long HASH_BASE = 0x2545F4914F6CDD1DL;
    /**
     * Below this many pairs in the sample, scaling it up says more about the sample than about the
     * population, and a caller should say so rather than report a number that looks as precise as any
     * other.
     */
    public static final long THIN_SAMPLE = 100_000;

    /**
     * One sketch per reader thread, so that nothing is locked on the reading path, held in a thread
     * local rather than looked up per pair: the pairs are counted in the billions, and a map keyed by
     * the thread id would cost a lookup and - the key being a boxed long - an allocation on every
     * single one of them. The sketches are collected as they are created so that they can be merged
     * afterwards, whichever threads happen to have run.
     */
    private final Collection<HLL> sketches = new ConcurrentLinkedQueue<>();
    private final ThreadLocal<HLL> sketch = ThreadLocal.withInitial(() -> {
        HLL created = newSketch();
        sketches.add(created);
        return created;
    });
    private final double sampleScale;
    private long sampledCount;

    /**
     * @param sampleScale what the count of the sample has to be multiplied by to stand for the whole;
     *                    one where the pass sees the entire population
     */
    public DistinctPairSketch(double sampleScale) {
        this.sampleScale = sampleScale;
    }

    /**
     * Records one pair. Safe to call from every reader thread at once.
     *
     * @param pairKey the pair, as {@link KMerIndexFilterHelper#combine} folds it into one word
     */
    public void record(long pairKey) {
        // Hashed here rather than taken as it comes: HyperLogLog reads the register index off the low
        // bits of its input and the leading zeros off the rest, so it needs input that is spread over
        // the whole word. KMerIndexFilterHelper.combine() folds the node index into the k-mer, but a
        // k-mer is two bits per base and leaves the top of the word untouched for every k below 32 -
        // registers would go unused and the count would come out wrong. The finalizer costs a few
        // multiplications against a pass over every reference sequence.
        sketch.get().addRaw(MurmurHash3DropIn.hash64(pairKey, HASH_BASE));
    }

    /**
     * Merges what the reader threads sketched and scales it back up to the population.
     * <p>
     * Call once the readers are done. It is not idempotent in the sense of being free: every call
     * merges again, so a caller wanting both this and {@link #getSampledCount()} should call this one
     * first.
     *
     * @return the estimated number of distinct pairs
     */
    public long estimate() {
        HLL merged = newSketch();
        for (HLL each : sketches) {
            merged.union(each);
        }
        sampledCount = merged.cardinality();
        return Math.round(sampledCount * sampleScale);
    }

    /**
     * @return how many distinct pairs the sample itself held, as of the last {@link #estimate()}
     */
    public long getSampledCount() {
        return sampledCount;
    }

    /** @return by what the sample was scaled up */
    public double getSampleScale() {
        return sampleScale;
    }

    /**
     * @return whether the sample was too thin for the estimate to be more than coarse
     */
    public boolean isThin() {
        return sampleScale > 1 && sampledCount < THIN_SAMPLE;
    }

    /** Releases the sketches. */
    public void clear() {
        sketches.clear();
    }

    /**
     * @return a fresh sketch, left to promote itself through the library's representations rather than
     *         forced to the fully materialised one: a population of few pairs is then counted exactly
     */
    private static HLL newSketch() {
        return new HLL(HLL_LOG2M, HLL_REGISTER_WIDTH);
    }
}
