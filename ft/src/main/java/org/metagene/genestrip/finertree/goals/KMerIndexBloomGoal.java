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
import org.metagene.genestrip.probfilter.BlockedBloomFilter;
import org.metagene.genestrip.probfilter.BloomFilter;
import org.metagene.genestrip.probfilter.ProbFilter;
import org.metagene.genestrip.probfilter.XORBloomFilter;
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
import org.metagene.genestrip.tax.TaxNodeSelection;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.function.LongToDoubleFunction;

/**
 * Builds the {@link ProbFilter} that records, for every relevant k-mer of the database, which leaf
 * tax node (genome) it originates from. The reference fastas are re-read and, for each k-mer that is
 * stored in the k-mer store below a taxon selected for refinement, the leaf tax node's store index is
 * registered in the filter. This index information is later used by the clustering phase to build a
 * refined (dendrogram) tax tree.
 * <p>
 * How large that filter has to be is the question {@link KMerIndexSizeGoal} answers by sketching the
 * pairs beforehand. That goal is a lazy dependency: it is an {@link ObjectGoal} and hence never made
 * on its own account, only when {@link ObjectGoal#get()} is called here - which happens exactly when
 * {@link FTConfigKey#FT_BLOOM_FILTER_SIZING} asks for it. Configured to the conservative bound, the goal
 * is not made at all and nothing is read twice.
 *
 * @param <P> the concrete {@link FTProject} type this goal operates on
 */
public class KMerIndexBloomGoal<P extends FTProject> extends AbstractKMerIndexGoal<ProbFilter, P> {
    private final ObjectGoal<Long, P> indexSizeGoal;

    /**
     * The false-positive rate beyond which a filter is refused rather than used. At the default FPP
     * the filter is the blocked one, built for one percent and measuring 1.3 % once it holds what it
     * was sized for; it reaches 4.4 % at one and a half times that and 9.4 % at twice. Refusing above
     * 5 % therefore tolerates an estimate that fell short by more than half - far beyond the fraction
     * of a percent HyperLogLog is out by - while keeping the answers the clustering reads from turning
     * into noise.
     */
    // The filter is refused when its false-positive rate has grown to this multiple of the rate it
    // would have at the count it was sized for. Relative and not an absolute rate, so that it keeps
    // its meaning whichever filter FT_INDEX_BLOOM_FILTER_FPP selects and whatever that filter spends
    // per key: the blocked one at the ten bits per key it uses today already answers wrongly 1.2 per
    // cent of the time when holding exactly what it was sized for, a plain one asked for 1e-4 answers
    // wrongly a hundredth as often, and an absolute bound would have to be re-derived for each.
    //
    // 1.75 is reached when the estimate fell about a fifth short. The estimate's own error is far
    // below that - the sketch is good to half a per cent and the sampling of kmerindexsize to about
    // the same - so percentage noise passes, while the failure this guards against is not a matter of
    // percentages: an estimate that misses the mark misses it by factors, as the conservative bound
    // exceeding the truth by orders of magnitude on a bacterial database shows.
    private static final double MAX_TRUSTED_FPP_FACTOR = 1.75;

    private ProbFilter filter;

    /**
     * Answers the false-positive rate {@link #filter} reaches at a given number of insertions, which
     * is what the trust check below is decided on.
     * <p>
     * Kept beside the filter, which is a {@link ProbFilter}, for the reason {@code FillBloomFilterGoal}
     * keeps its own: only a bloom filter can answer it. Here it is a function rather than a second
     * typed field because the two kinds {@link #createFilter(long)} chooses between - blocked and
     * plain - declare {@code getFpp(long)} separately and share no supertype that has it.
     */
    private LongToDoubleFunction fppAt;

    /**
     * Creates the goal that builds the k-mer index bloom filter.
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
     * @param indexSizeGoal    goal estimating the number of filter entries, consulted only when the
     *                         configuration asks for its estimate
     * @param deps             further goals this goal depends on
     */
    @SafeVarargs
    public KMerIndexBloomGoal(P project, ExecutionContext bundle, ObjectGoal<Set<RefSeqCategory>, P> categoriesGoal,
                              ObjectGoal<TaxNodeSelection, P> taxNodesGoal,
                              ObjectGoal<TaxTree, P> taxTreeGoal, RefSeqFnaFilesDownloadGoal fnaFilesGoal,
                              ObjectGoal<Map<File, TaxTree.TaxIdNode>, P> additionalGoal,
                              ObjectGoal<AccessionMap, P> accessionMapGoal, ObjectGoal<Database, P> storeGoal,
                              ObjectGoal<Long, P> indexSizeGoal, Goal<P>... deps) {
        super(project, FTGoalKey.KMER_INDEX_BLOOM, bundle, categoriesGoal, taxNodesGoal, taxTreeGoal, fnaFilesGoal,
                additionalGoal, accessionMapGoal, storeGoal, Goal.append(deps, indexSizeGoal));
        this.indexSizeGoal = indexSizeGoal;
    }

    @Override
    protected boolean record(long hash) {
        // Lock-free combined membership-check-and-insert: the filter sets its bits atomically, so
        // concurrent readers never need a lock on it.
        return filter.putLong(hash);
    }

    /**
     * Determines the relevant tax nodes for refinement, sizes and allocates the {@link ProbFilter},
     * and re-reads the fastas to populate it with per-k-mer leaf-node indices. The completed filter is
     * published as this goal's result.
     */
    @Override
    protected void doMakeThis() {
        GSConfigKey.BloomFilterSizing sizing =
                (GSConfigKey.BloomFilterSizing) configValue(FTConfigKey.FT_BLOOM_FILTER_SIZING);
        try {
            prepare();
            long conservative = conservativeEstimate();
            if (getLogger().isInfoEnabled()) {
                getLogger().info("Maximum expected filter entries: " + conservative);
            }
            long sizingCount = conservative;
            if (sizing != GSConfigKey.BloomFilterSizing.UPPER_BOUND) {
                // Here and only here the sizing goal is asked for its value, and asking is what makes
                // it run. The extra pass over the sequences it costs is the price of a filter sized
                // for what the index really holds rather than for what it could conceivably hold.
                long estimated = indexSizeGoal.get();
                if (getLogger().isInfoEnabled()) {
                    getLogger().info("Sizing the k-mer index filter for " + estimated + " entries, the"
                            + " estimated number of distinct ones; the conservative bound is " + conservative
                            + ", larger by a factor of "
                            + (estimated > 0 ? ((double) conservative) / estimated : Double.NaN));
                }
                sizingCount = estimated;
            } else if (getLogger().isInfoEnabled()) {
                getLogger().info("Sizing the k-mer index filter for " + conservative
                        + " entries, the conservative bound, as the configuration asks for.");
            }
            long entries = onePass(sizingCount);
            double reachedFpp = fppAt.applyAsDouble(entries);
            double designFpp = fppAt.applyAsDouble(sizingCount);
            double maxTrustedFpp = designFpp * MAX_TRUSTED_FPP_FACTOR;
            if (reachedFpp > maxTrustedFpp) {
                // What matters is not whether the filter held more than it was built for - a single
                // k-mer over would be neither here nor there - but whether its false-positive rate has
                // grown enough to corrupt what is read from it: every false positive is a k-mer the
                // clustering takes for shared between two species that do not share it.
                if (sizing == GSConfigKey.BloomFilterSizing.AUTO) {
                    if (getLogger().isWarnEnabled()) {
                        getLogger().warn("The k-mer index filter was sized for " + sizingCount
                                + " entries but received " + entries + ", bringing its false-positive rate to "
                                + String.format("%.4f", reachedFpp) + " where anything above "
                                + String.format("%.4f", maxTrustedFpp) + " is refused - "
                                + MAX_TRUSTED_FPP_FACTOR + " times the " + String.format("%.4f", designFpp)
                                + " it was built for -, so the estimate fell short."
                                + " Reading the sequences again with the conservative bound of "
                                + conservative + ", as '" + FTConfigKey.FT_BLOOM_FILTER_SIZING.getName() + "="
                                + GSConfigKey.BloomFilterSizing.AUTO.getName() + "' asks for. Set it to '"
                                + GSConfigKey.BloomFilterSizing.UPPER_BOUND.getName()
                                + "' to go straight there next time.");
                    }
                    releaseAfterPass();
                    prepare();
                    entries = onePass(conservative);
                } else {
                    throw new IllegalStateException("The k-mer index filter was sized for " + sizingCount
                            + " entries but received " + entries + ", bringing its false-positive rate to "
                            + String.format("%.4f", reachedFpp) + " where anything above "
                            + String.format("%.4f", maxTrustedFpp) + " is refused - " + MAX_TRUSTED_FPP_FACTOR
                            + " times the " + String.format("%.4f", designFpp) + " it was built for -."
                            + " The clustering reading such a filter would take the surplus"
                            + " for shared k-mers. Set '" + FTConfigKey.FT_BLOOM_FILTER_SIZING.getName()
                            + "=" + GSConfigKey.BloomFilterSizing.UPPER_BOUND.getName() + "' (or '"
                            + GSConfigKey.BloomFilterSizing.AUTO.getName() + "', which falls back to it by"
                            + " itself) and run the goal again; the conservative bound of " + conservative
                            + " cannot be exceeded.");
                }
            }
            set(filter);
            if (getLogger().isInfoEnabled()) {
                getLogger().info("Actual entries: " + entries);
            }
        } catch (IOException e) {
            throw new RuntimeException(e);
        } finally {
            releaseAfterPass();
        }
    }

    /**
     * Returns the conservative bound on the number of (k-mer, child index) pairs: the number reached
     * if every k-mer stored on a node occurred in the genome of every one of its direct subnodes,
     * plus one for the OTHER slot.
     * <p>
     * Direct subnodes and not the whole subtree, because that is where the pairs are recorded: the
     * index pass maps a genome's leaf onto the direct child of the refined node it lies under (see
     * {@code AbstractKMerIndexGoal#childIndexUnder}), so a node cannot produce more distinct pairs
     * than it has children, plus OTHER. The count below is therefore a bound in fact, not merely in
     * intent - which it was not while pairs went in under the leaf itself, deep in the tree, and a
     * node with one child could produce thousands.
     * <p>
     * The bound is exact only in the sense of never being exceeded. How far above the truth it lies
     * depends on how much the genomes of one taxon really share, which is what {@code kmerindexsize}
     * exists to find out. Reading the counts again costs nothing: the store caches them.
     *
     * @return the conservative bound on the number of pairs
     * @throws IllegalStateException if {@link #prepare()} has not run
     */
    protected long conservativeEstimate() {
        if (kmerStore == null) {
            throw new IllegalStateException("conservativeEstimate() needs the store that prepare() loads.");
        }
        long[] counter = new long[1];
        kmerStore.getNKmersPerTaxid().forEach((s, aLong) -> {
            if (s != null) {
                if (isRefinementNode(s)) {
                    // Conservative estimate: k-mer could be in genome of every subnode, i.e. species...
                    // "+ 1" is for nodes not included in the database but below a rank to refine.
                    counter[0] += aLong * (s.getNumberOfSubNodes() + 1);
                }
            }
        });
        return counter[0];
    }

    /**
     * Allocates a filter of the given size and reads all sequences once into it.
     *
     * @param sizingCount the number of entries to size the filter for
     * @return the number of entries actually recorded
     * @throws IOException if a sequence file cannot be read
     */
    protected long onePass(long sizingCount) throws IOException {
        createFilter(sizingCount);
        if (getLogger().isInfoEnabled()) {
            getLogger().info("Filter size in MB: " + (filter.getBitSize() / 8 / 1024 / 1024));
        }
        readFastas();
        long entries = 0;
        for (MyFastaReader reader : readers) {
            // The readers are done, so the k-mers still buffered from their last (partial) batch are
            // recorded here, single-threaded, before the filter is published and the entries counted.
            reader.flushBatch();
            entries += reader.getEntries();
        }
        return entries;
    }

    /**
     * Allocates the filter for {@link FTConfigKey#FT_INDEX_BLOOM_FILTER_FPP} and sets {@link #fppAt}
     * to answer for it, choosing the kind the way {@code BloomIndexGoal} and
     * {@code AbstractKMerStore.createOptimizedFilter} choose theirs: the blocked filter at the default
     * rate, a plain one below it.
     * <p>
     * The blocked filter is the faster of the two - about a quarter of the time per lookup once both
     * are past the last-level cache - but it sets a fixed four bits per key, which is near the optimum
     * at the ten bits per key it is given and far from it below: reaching 4.5e-4 costs it 28 bits per
     * entry where the optimal sizing needs 17. Since the whole point of asking for a lower rate here is
     * the {@code coversAll} placement in {@code ftupdatedb}, which one false positive out of a node's
     * C child lookups defeats, the memory is the side worth being efficient on.
     * <p>
     * Both kinds fill safely from the reader threads - the blocked one locks the bucket owning the
     * word, {@code LargeBitVector} the bucket behind the bit - which {@link #record(long)} relies on.
     *
     * @param sizingCount the number of entries to size the filter for
     */
    protected void createFilter(long sizingCount) {
        double fpp = doubleConfigValue(FTConfigKey.FT_INDEX_BLOOM_FILTER_FPP);
        if (fpp >= BlockedBloomFilter.DEFAULT_FPP) {
            BlockedBloomFilter blocked = new BlockedBloomFilter(sizingCount);
            filter = blocked;
            fppAt = blocked::getFpp;
        } else {
            BloomFilter bloom = booleanConfigValue(GSConfigKey.XOR_BLOOM_HASH) ?
                    new XORBloomFilter(fpp, sizingCount) : new BloomFilter(fpp, sizingCount);
            filter = bloom;
            fppAt = bloom::getFpp;
        }
        if (getLogger().isInfoEnabled()) {
            getLogger().info("K-mer index filter: " + filter.getClass().getSimpleName() + " at an FPP of "
                    + fppAt.applyAsDouble(sizingCount) + " for " + sizingCount + " entries.");
        }
    }
}
