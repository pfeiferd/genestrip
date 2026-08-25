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
import org.metagene.genestrip.finertree.FTGoalKey;
import org.metagene.genestrip.finertree.FTProject;
import org.metagene.genestrip.finertree.probfilter.DistinctPairSketch;
import org.metagene.genestrip.finertree.probfilter.KMerIndexFilterHelper;
import org.metagene.genestrip.goals.refseq.RefSeqFnaFilesDownloadGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.refseq.AbstractRefSeqFastaReader;
import org.metagene.genestrip.refseq.AccessionMap;
import org.metagene.genestrip.refseq.RefSeqCategory;
import org.metagene.genestrip.store.Database;
import org.metagene.genestrip.tax.SmallTaxTree;
import org.metagene.genestrip.tax.TaxTree;
import org.metagene.genestrip.tax.TaxNodeSelection;

import java.io.File;
import java.io.IOException;
import java.util.Map;
import java.util.Set;

/**
 * Estimates how many distinct (*k*-mer, leaf) pairs {@link DBQualityCountsGoal} will produce, by
 * reading the sequences once and sketching the pairs with HyperLogLog.
 * <p>
 * It exists because the conservative bound cannot size that goal's filter for every database. The
 * bound sums, for every leaf, the *k*-mers stored along its path to the root, i.e. it assumes each of
 * them to occur in every leaf below its node. That holds for a database whose *k*-mers sit close to
 * the leaves and fails badly for one whose weight is at a single high node -- `strepto' keeps 222
 * million *k*-mers at the genus, on the path of every one of its several hundred leaves, and the bound
 * comes to some 10^11 entries. At ten bits each the filter alone would want tens of gigabytes, and the
 * run died in {@code newLargeGrid()} before it read a single base.
 * <p>
 * The same trade {@code kmerindexsize} makes for the index filter, through the same
 * {@link DistinctPairSketch}: one more pass over the sequences for a filter of the right size.
 *
 * @param <P> the concrete FT project type
 */
public class DBQualitySizeGoal<P extends FTProject> extends AbstractDBQualityGoal<Long, P> {
    /**
     * What the count of the sample has to be multiplied by to stand for the whole -- one, because this
     * pass looks at every k-mer the database holds and nothing has to be scaled back up.
     * {@code kmerindexsize} samples one k-mer in {@code ftKMerIndexSizeSampling} and scales by that;
     * the constant is named here so that the two goals read alike and the difference between them is
     * visible rather than implied.
     */
    private static final int SAMPLE_SCALE = 1;

    /** Counts the distinct pairs of this pass; held for its duration. */
    private DistinctPairSketch sketch;

    /**
     * Creates the goal.
     *
     * @param project          the FT project
     * @param key              the goal key
     * @param bundle           the execution context used to run the fasta readers
     * @param categoriesGoal   the goal providing the RefSeq categories to read
     * @param taxNodesGoal     the goal providing the tax nodes to be included
     * @param fnaFilesGoal     the goal providing the downloaded genomic fasta files
     * @param additionalGoal   the goal providing additional fasta files mapped to tax nodes
     * @param accessionMapGoal the goal providing the accession-to-tax-node map
     * @param storeGoal        the goal providing the loaded database
     * @param deps             further goals this goal depends on
     */
    @SafeVarargs
    public DBQualitySizeGoal(P project, FTGoalKey key, ExecutionContext bundle, ObjectGoal<Set<RefSeqCategory>, P> categoriesGoal,
                             ObjectGoal<TaxNodeSelection, P> taxNodesGoal,
                             RefSeqFnaFilesDownloadGoal fnaFilesGoal,
                             ObjectGoal<Map<File, TaxTree.TaxIdNode>, P> additionalGoal,
                             ObjectGoal<AccessionMap, P> accessionMapGoal, ObjectGoal<Database, P> storeGoal,
                             Goal<P>... deps) {
        super(project, key, bundle, categoriesGoal, taxNodesGoal, fnaFilesGoal, additionalGoal,
                accessionMapGoal, storeGoal, deps);
    }

    /**
     * Reads the sequences once and stores the number of distinct pairs the sketch counted.
     */
    @Override
    protected void doMakeThis() {
        try {
            // Created before anything that can fail, so that the finally block below always has a
            // sketch to clear.
            sketch = new DistinctPairSketch(SAMPLE_SCALE);
            // The store is what this pass needs from prepare(): a k-mer the database does not hold
            // forms no pair.
            prepare();
            readFastas();
            for (MyFastaReader reader : readers) {
                // The readers are done, so the pairs still buffered from their last (partial) batch are
                // recorded here, single-threaded, before the sketches are merged.
                reader.flushBatch();
            }
            long estimated = sketch.estimate();
            set(estimated);
            if (getLogger().isInfoEnabled()) {
                // How this compares to the conservative bound is logged by the goal that uses both,
                // which is where the two numbers meet anyway.
                getLogger().info("Estimated distinct filter entries: " + estimated);
            }
            // Neither the sample suffix nor the thin-sample warning of kmerindexsize appears here:
            // SAMPLE_SCALE is 1, so the count is the population and not an estimate scaled up from a
            // part of it.
        } catch (IOException e) {
            throw new RuntimeException(e);
        } finally {
            sketch.clear();
            releaseAfterPass();
        }
    }

    @Override
    protected MyFastaReader newReader(AbstractRefSeqFastaReader.StringLong2DigitTrie regionsPerTaxid) {
        return new SketchingReader(regionsPerTaxid);
    }

    /** Reader of this pass: it sketches every pair and keeps no tallies, there being none to keep. */
    private class SketchingReader extends MyFastaReader {
        /**
         * Creates the reader.
         *
         * @param regionsPerTaxid the trie counting regions per tax id
         */
        SketchingReader(AbstractRefSeqFastaReader.StringLong2DigitTrie regionsPerTaxid) {
            super(regionsPerTaxid);
        }

        /**
         * Sketches one pair. How many distinct ones there are is the whole question of this pass, so
         * nothing else happens to it, and there is no filter yet to deduplicate against -- its size is
         * what this pass is for.
         *
         * @param kmer       the k-mer
         * @param leafPos    the position of the leaf it was read in
         * @param storedNode the node the database stores it at, which this pass does not look at
         * @return always {@code true}; the sketch cannot say whether a pair was new, and the count
         * that matters is the one it gives at the end
         */
        @Override
        protected boolean count(long kmer, int leafPos, SmallTaxTree.SmallTaxIdNode storedNode) {
            sketch.record(KMerIndexFilterHelper.combine(kmer, leafPos));
            entries++;
            return true;
        }
    }
}
