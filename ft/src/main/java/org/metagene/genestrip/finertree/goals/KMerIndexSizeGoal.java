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
package org.metagene.genestrip.finertree.goals;

import org.metagene.genestrip.ExecutionContext;
import org.metagene.genestrip.GSConfigKey;
import org.metagene.genestrip.finertree.FTConfigKey;
import org.metagene.genestrip.finertree.FTGoalKey;
import org.metagene.genestrip.finertree.FTProject;
import org.metagene.genestrip.goals.refseq.RefSeqFnaFilesDownloadGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.refseq.AccessionMap;
import org.metagene.genestrip.refseq.RefSeqCategory;
import org.metagene.genestrip.store.Database;
import org.metagene.genestrip.tax.TaxTree;
import org.metagene.genestrip.util.KMerSampling;
import org.metagene.genestrip.util.MurmurHash3DropIn;

import net.agkn.hll.HLL;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.Collection;
import java.util.concurrent.ConcurrentLinkedQueue;

/**
 * Estimates how many (k-mer, leaf index) pairs the k-mer index filter will have to hold, by reading
 * the reference sequences once and sketching the pairs with HyperLogLog instead of storing them.
 * <p>
 * The alternative is the conservative bound that {@code kmerindexbloom} computes for itself, which
 * assumes every k-mer of a taxon to occur in the genome of every one of its subnodes. Genomes of one
 * taxon share most of their k-mers, so that bound can exceed the truth by orders of magnitude, and a
 * filter built to it is larger by the same factor. The price of knowing better is this goal: one
 * more pass over every sequence.
 * <p>
 * The sketch is kept per reader thread and merged at the end, so nothing is locked while reading -
 * see {@link org.metagene.genestrip.goals.refseq.FillSizeGoal} which counts the database's k-mers
 * the same way.
 * <p>
 * By default only one k-mer in
 * {@link org.metagene.genestrip.finertree.FTConfigKey#FT_KMER_INDEX_SIZE_SAMPLING} is looked at, which
 * is what makes the pass affordable: the ones left out are skipped before the store is consulted. The
 * sample is drawn by the hash of the k-mer and not by position - see {@link #considerKMer(long)} for
 * why the obvious rule is the wrong one.
 *
 * @param <P> the concrete {@link FTProject} type this goal operates on
 */
public class KMerIndexSizeGoal<P extends FTProject> extends AbstractKMerIndexGoal<Long, P> {
    // Sizing of the sketches: 2^15 registers of 6 bits for a relative error of about half a percent,
    // and one hash base for all readers, since sketches can only be merged when they hashed alike.
    private static final int HLL_LOG2M = 15;
    private static final int HLL_REGISTER_WIDTH = 6;

    // The seed all sketches hash with. Merging two sketches is only meaningful when they hashed
    // their input alike, so this is a constant and not drawn per reader.
    private static final long HASH_BASE = 0x2545F4914F6CDD1DL;

    // Below this many pairs in the sample, scaling it up says more about the sample than about the
    // index, and the run says so rather than reporting a number that looks as precise as any other.
    private static final long THIN_SAMPLE = 100_000;

    // One sketch per reader thread, so that nothing is locked on the reading path, held in a thread
    // local rather than looked up per pair: the pairs are counted in the billions, and a map keyed
    // by the thread id would cost a lookup and - the key being a boxed long - an allocation on every
    // single one of them. The sketches are collected as they are created so that they can be merged
    // afterwards, whichever threads happen to have run.
    private final Collection<HLL> sketches = new ConcurrentLinkedQueue<>();

    // The bound the sampling decision compares against, derived from the configured rate so that no
    // division has to happen per k-mer. Final: the question is put for every k-mer of every reference
    // sequence, from every reader thread.
    private final long samplingThreshold;

    // What the count of the sample has to be multiplied by to stand for the whole.
    //
    // Not simply the configured rate, because the k-mers reaching this goal have already been sampled
    // once: the database keeps one k-mer in kMerSampling, by the very same rule, and this pass reads the
    // sequences the same way, so its population is that sample and not all k-mers there are. Sampling
    // a sample with the same rule nests rather than multiplies - the finer threshold selects a subset
    // of the coarser one - so of a population reduced to one in s, a threshold for one in q keeps a
    // fraction of s/max(s, q). That is exact and needs no assumption about two hash functions being
    // independent. With the default step size of 1 it comes to the configured rate, as one would
    // expect; with a rate no finer than the step size it comes to 1, the sample then being the whole
    // population.
    private final double sampleScale;

    private final ThreadLocal<HLL> sketch = ThreadLocal.withInitial(() -> {
        HLL created = newSketch();
        sketches.add(created);
        return created;
    });

    /**
     * Creates the goal.
     *
     * @param project          the project this goal belongs to
     * @param bundle           the execution context providing the number of worker threads
     * @param categoriesGoal   goal supplying the RefSeq categories to process
     * @param taxNodesGoal     goal supplying the tax nodes to be considered
     * @param taxTreeGoal      goal supplying the (large) tax tree
     * @param fnaFilesGoal     goal supplying the downloaded RefSeq fna files
     * @param additionalGoal   goal supplying additional fasta files mapped to their tax node
     * @param accessionMapGoal goal supplying the accession-to-tax-node mapping
     * @param storeGoal        goal supplying the database (and its k-mer store) to index
     * @param deps             further goals this goal depends on
     */
    @SafeVarargs
    public KMerIndexSizeGoal(P project, ExecutionContext bundle, ObjectGoal<Set<RefSeqCategory>, P> categoriesGoal,
                             ObjectGoal<Set<TaxTree.TaxIdNode>, P> taxNodesGoal,
                             ObjectGoal<TaxTree, P> taxTreeGoal, RefSeqFnaFilesDownloadGoal fnaFilesGoal,
                             ObjectGoal<Map<File, TaxTree.TaxIdNode>, P> additionalGoal,
                             ObjectGoal<AccessionMap, P> accessionMapGoal, ObjectGoal<Database, P> storeGoal,
                             Goal<P>... deps) {
        super(project, FTGoalKey.KMER_INDEX_SIZE, bundle, categoriesGoal, taxNodesGoal, taxTreeGoal, fnaFilesGoal,
                additionalGoal, accessionMapGoal, storeGoal, deps);
        int sampling = intConfigValue(FTConfigKey.FT_KMER_INDEX_SIZE_SAMPLING);
        int kMerSampling = intConfigValue(GSConfigKey.KMER_SAMPLING);
        samplingThreshold = KMerSampling.thresholdForOneIn(sampling);
        sampleScale = ((double) Math.max(kMerSampling, sampling)) / kMerSampling;
    }

    /**
     * Creates a sketch of the sizing all readers share, so that theirs can be merged.
     *
     * @return a new, empty sketch
     */
    protected static HLL newSketch() {
        // Left to promote itself through the library's representations rather than forced to the fully
        // materialised one: an index of few pairs is then counted exactly.
        return new HLL(HLL_LOG2M, HLL_REGISTER_WIDTH);
    }

    @Override
    protected boolean considerKMer(long kmer) {
        return KMerSampling.isSampled(kmer, samplingThreshold);
    }

    @Override
    protected boolean record(long hash) {
        // Hashed here rather than taken as it comes: HyperLogLog reads the register index off the low
        // bits of its input and the leading zeros off the rest, so it needs input that is spread over
        // the whole word. KMerIndexFilterHelper.combine() folds the leaf index into the k-mer, but a
        // k-mer is two bits per base and leaves the top of the word untouched for every k below 32 -
        // registers would go unused and the count would come out wrong. The finalizer costs a few
        // multiplications against a pass over every reference sequence.
        sketch.get().addRaw(MurmurHash3DropIn.hash64(hash, HASH_BASE));
        // The sketch cannot say whether this pair was new, and nothing here needs to know: the count
        // that matters is the one it gives at the end.
        return true;
    }

    @Override
    protected void doMakeThis() {
        try {
            prepare();
            readFastas();
            for (MyFastaReader reader : readers) {
                // The readers are done, so the pairs still buffered from their last (partial) batch are
                // recorded here, single-threaded, before the sketches are merged.
                reader.flushBatch();
            }
            HLL merged = newSketch();
            for (HLL each : sketches) {
                merged.union(each);
            }
            long inSample = merged.cardinality();
            long estimated = Math.round(inSample * sampleScale);
            set(estimated);
            if (getLogger().isInfoEnabled()) {
                // How this compares to the conservative bound is logged by the goal that uses both,
                // which is where the two numbers meet anyway.
                getLogger().info("Estimated distinct filter entries: " + estimated
                        + (sampleScale == 1 ? "" : " (from " + inSample + " counted in a sample of one k-mer in "
                        + String.format("%.0f", sampleScale) + " of those the database keeps)"));
            }
            if (sampleScale > 1 && inSample < THIN_SAMPLE && getLogger().isWarnEnabled()) {
                getLogger().warn("Only " + inSample + " pair(s) were sampled, so this estimate is a coarse one."
                        + " Set " + FTConfigKey.FT_KMER_INDEX_SIZE_SAMPLING.getName() + " to 1 to count them all;"
                        + " for an index this small, reading everything costs little.");
            }
        } catch (IOException e) {
            throw new RuntimeException(e);
        } finally {
            sketches.clear();
            releaseAfterPass();
        }
    }
}
