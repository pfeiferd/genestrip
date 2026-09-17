package org.metagene.genestrip.finertree.probfilter;

import org.junit.Test;

import java.util.Random;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

/**
 * Pins what {@code dbqualcounts} and {@code kmerindexsize} both rely on: that the sketch counts the
 * distinct pairs it is shown, closely enough to size a filter by, whichever thread showed them.
 * <p>
 * The estimate decides how large a Bloom filter is allocated. Too small and the filter's false
 * positives silently drop counts; too large and the run dies before it reads a base, which is what
 * the conservative bound did on the paper's \texttt{strepto} database. Neither failure announces
 * itself, so the properties are pinned here rather than trusted.
 */
public class DistinctPairSketchTest {
    /** HyperLogLog at 2^15 registers is specified for about half a percent; allow four times that. */
    private static final double TOLERANCE = 0.02;

    @Test
    public void testCountsDistinctPairsOfAFullPass() {
        DistinctPairSketch sketch = new DistinctPairSketch(1);
        int distinct = 500_000;
        for (int i = 0; i < distinct; i++) {
            sketch.record(KMerIndexFilterHelper.combine(i, i % 997));
        }
        long estimate = sketch.estimate();
        assertEquals("the estimate is what sizes the filter", distinct, estimate,
                distinct * TOLERANCE);
    }

    @Test
    public void testRepeatedPairsAreNotCountedTwice() {
        DistinctPairSketch sketch = new DistinctPairSketch(1);
        int distinct = 200_000;
        // Every pair five times over, as a k-mer shared by five genomes of one leaf would arrive.
        for (int round = 0; round < 5; round++) {
            for (int i = 0; i < distinct; i++) {
                sketch.record(KMerIndexFilterHelper.combine(i, 3));
            }
        }
        assertEquals("a pair seen five times is one pair", distinct, sketch.estimate(),
                distinct * TOLERANCE);
    }

    @Test
    public void testTheSameKMerUnderDifferentLeavesCountsSeparately() {
        // The whole point of the pair: one k-mer carried by three leaves is three entries in the
        // filter, and a sizing that collapsed them would undersize it by the branching degree.
        DistinctPairSketch sketch = new DistinctPairSketch(1);
        int kmers = 100_000;
        for (int i = 0; i < kmers; i++) {
            sketch.record(KMerIndexFilterHelper.combine(i, 1));
            sketch.record(KMerIndexFilterHelper.combine(i, 2));
            sketch.record(KMerIndexFilterHelper.combine(i, 3));
        }
        assertEquals(3L * kmers, sketch.estimate(), 3.0 * kmers * TOLERANCE);
    }

    @Test
    public void testThreadsAgreeWithOneThread() throws Exception {
        // The readers sketch on their own threads and the sketches are merged at the end. A merge
        // that lost a thread's share would undersize the filter in exact proportion to the threads,
        // which is the kind of error that looks like a tuning problem rather than a defect.
        int distinct = 400_000;
        DistinctPairSketch single = new DistinctPairSketch(1);
        for (int i = 0; i < distinct; i++) {
            single.record(KMerIndexFilterHelper.combine(i, i % 13));
        }
        long expected = single.estimate();

        DistinctPairSketch shared = new DistinctPairSketch(1);
        int threads = 8;
        ExecutorService pool = Executors.newFixedThreadPool(threads);
        CountDownLatch done = new CountDownLatch(threads);
        for (int t = 0; t < threads; t++) {
            final int offset = t;
            pool.execute(() -> {
                try {
                    for (int i = offset; i < distinct; i += threads) {
                        shared.record(KMerIndexFilterHelper.combine(i, i % 13));
                    }
                } finally {
                    done.countDown();
                }
            });
        }
        assertTrue("the sketching threads did not finish", done.await(60, TimeUnit.SECONDS));
        pool.shutdown();
        assertEquals("merging the per-thread sketches must give what one thread gives",
                expected, shared.estimate(), expected * TOLERANCE);
    }

    @Test
    public void testASampleIsScaledBackUp() {
        // kmerindexsize looks at one k-mer in n and scales the count up. The scale is the sketch's
        // business, and getting it wrong scales the filter with it.
        DistinctPairSketch sketch = new DistinctPairSketch(16);
        int inSample = 300_000;
        for (int i = 0; i < inSample; i++) {
            sketch.record(KMerIndexFilterHelper.combine(i, 7));
        }
        assertEquals(16L * inSample, sketch.estimate(), 16.0 * inSample * TOLERANCE);
        assertEquals(inSample, sketch.getSampledCount(), inSample * TOLERANCE);
        assertFalse("a sample of 300,000 pairs is not thin", sketch.isThin());
    }

    @Test
    public void testAThinSampleSaysSo() {
        DistinctPairSketch sketch = new DistinctPairSketch(64);
        for (int i = 0; i < 1000; i++) {
            sketch.record(KMerIndexFilterHelper.combine(i, 1));
        }
        sketch.estimate();
        assertTrue("1,000 pairs scaled by 64 is a coarse estimate and should be reported as one",
                sketch.isThin());
    }

    @Test
    public void testAnUnsampledPassIsNeverThin() {
        DistinctPairSketch sketch = new DistinctPairSketch(1);
        sketch.record(KMerIndexFilterHelper.combine(1, 1));
        sketch.estimate();
        assertFalse("a pass that saw everything is exact, however few pairs there were",
                sketch.isThin());
    }

    @Test
    public void testEmptyPassCountsNothing() {
        DistinctPairSketch sketch = new DistinctPairSketch(1);
        assertEquals(0, sketch.estimate());
    }

    @Test
    public void testRandomPairsAreCountedAsWell() {
        // The pairs of a real pass are neither consecutive nor evenly spread, and HyperLogLog reads
        // its register index off the low bits: a key whose entropy sits elsewhere would use a
        // fraction of the registers. DistinctPairSketch.record() hashes for that reason, and this
        // shows the hashing does its job on keys that are not tidy.
        DistinctPairSketch sketch = new DistinctPairSketch(1);
        Random random = new Random(4711);
        int distinct = 300_000;
        for (int i = 0; i < distinct; i++) {
            sketch.record(KMerIndexFilterHelper.combine(random.nextLong() >>> 2, random.nextInt(500)));
        }
        assertEquals(distinct, sketch.estimate(), distinct * TOLERANCE);
    }
}
