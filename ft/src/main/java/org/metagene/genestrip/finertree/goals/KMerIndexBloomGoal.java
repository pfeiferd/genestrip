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

import it.unimi.dsi.fastutil.objects.Object2LongMap;
import org.metagene.genestrip.ExecutionContext;
import org.metagene.genestrip.GSConfigKey;
import org.metagene.genestrip.finertree.FTConfigKey;
import org.metagene.genestrip.finertree.FTGoalKey;
import org.metagene.genestrip.finertree.FTProject;
import org.metagene.genestrip.finertree.bloom.XORKMerIndexBloomFilter;
import org.metagene.genestrip.goals.refseq.FastaReaderGoal;
import org.metagene.genestrip.goals.refseq.RefSeqFnaFilesDownloadGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.refseq.AbstractRefSeqFastaReader;
import org.metagene.genestrip.refseq.AbstractStoreFastaReader;
import org.metagene.genestrip.refseq.AccessionMap;
import org.metagene.genestrip.refseq.RefSeqCategory;
import org.metagene.genestrip.store.Database;
import org.metagene.genestrip.store.KMerSortedArray;
import org.metagene.genestrip.tax.Rank;
import org.metagene.genestrip.tax.SmallTaxTree;
import org.metagene.genestrip.tax.TaxIdCollector;
import org.metagene.genestrip.tax.TaxTree;

import java.io.File;
import java.io.IOException;
import java.util.*;

public class KMerIndexBloomGoal<P extends FTProject> extends FastaReaderGoal<XORKMerIndexBloomFilter, P> implements Goal.LogHeapInfo {
    public static final int OTHER_VALUE = Integer.MAX_VALUE;

    private final ObjectGoal<AccessionMap, P> accessionMapGoal;
    private final ObjectGoal<Database, P> storeGoal;
    private final ObjectGoal<TaxTree, P> taxTreeGoal;
    private final boolean multiThreading;
    /*
    private final boolean[] ranksToRefine;
    private final List<String> taxidsToRefine;
     */
    protected final List<FTConfigKey.RefinementPosition> refinementIntervals;

    private KMerSortedArray<SmallTaxTree.SmallTaxIdNode> kMerSortedArray;
    private SmallTaxTree smallTaxTree;
    private XORKMerIndexBloomFilter filter;
    private Set<TaxTree.TaxIdNode> relevantNodes;

    @SafeVarargs
    public KMerIndexBloomGoal(P project, ExecutionContext bundle, ObjectGoal<Set<RefSeqCategory>, P> categoriesGoal,
                              ObjectGoal<Set<TaxTree.TaxIdNode>, P> taxNodesGoal,
                              ObjectGoal<TaxTree, P> taxTreeGoal, RefSeqFnaFilesDownloadGoal fnaFilesGoal,
                              ObjectGoal<Map<File, TaxTree.TaxIdNode>, P> additionalGoal,
                              ObjectGoal<AccessionMap, P> accessionMapGoal, ObjectGoal<Database, P> storeGoal,
                              Goal<P>... deps) {
        super(project, FTGoalKey.KMER_INDEX_BLOOM, bundle, categoriesGoal, taxNodesGoal, fnaFilesGoal, additionalGoal, Goal.append(deps, taxTreeGoal, accessionMapGoal, storeGoal));
        this.storeGoal = storeGoal;
        this.accessionMapGoal = accessionMapGoal;
        this.taxTreeGoal = taxTreeGoal;
        multiThreading = bundle.getThreads() > 0;
        refinementIntervals = (List<FTConfigKey.RefinementPosition>) configValue(FTConfigKey.REFINEMENT_POSITIONS);
    }

    @Override
    protected void doMakeThis() {
        try {
            TaxIdCollector collector = new TaxIdCollector(taxTreeGoal.get());
            // Quite inefficient but should be good enough at this place.
            Set<TaxTree.TaxIdNode> nodesWithRank = new HashSet<>();
            for (TaxTree.TaxIdNode node : taxNodesGoal.get()) {
                while (node != null) {
                    if (FTConfigKey.RefinementPosition.getMatchingNodeFor(node, refinementIntervals) != null) {
                        nodesWithRank.add(node);
                    }
                    node = node.getParent();
                }
            }
            // Include subnodes from all nodes where we have ranks to refine:
            relevantNodes = collector.withDescendants(nodesWithRank, (Rank) configValue(GSConfigKey.RANK_COMPLETION_DEPTH));
            smallTaxTree = storeGoal.get().getTaxTree();
            kMerSortedArray = storeGoal.get().convertKMerStore();
            Object2LongMap<SmallTaxTree.SmallTaxIdNode> stats = kMerSortedArray.getNKmersPerTaxid();
            long[] counter = new long[1];
            stats.forEach((s, aLong) -> {
                if (s != null && s.getSubNodes() != null) {
                    if (FTConfigKey.RefinementPosition.getMatchingNodeFor(s, refinementIntervals) != null) {
                        // Conservative estimate: k-mer could be in genome of every subnode, i.e. species...
                        // "+ 1" is for nodes not included in the database but below a rank to refine.
                        counter[0] += aLong * (s.getSubNodes().length + 1);
                    }
                }
            });
            if (getLogger().isInfoEnabled()) {
                getLogger().info("Maximum expected filter entries: " + counter[0]);
            }
            filter = new XORKMerIndexBloomFilter(doubleConfigValue(FTConfigKey.FT_BLOOM_FILTER_FPP));
            long bitSize = filter.ensureExpectedSize(counter[0], false);
            if (getLogger().isInfoEnabled()) {
                getLogger().info("Filter size in MB: " + (bitSize / 8 / 1024 / 1024));
            }
            readFastas();
            set(filter);
            if (getLogger().isInfoEnabled()) {
                getLogger().info("Actual filter entries: " + filter.getEntries());
            }
        } catch (IOException e) {
            throw new RuntimeException(e);
        } finally {
            cleanUpThreads();
        }
    }

    @Override
    protected AbstractStoreFastaReader createFastaReader(AbstractRefSeqFastaReader.StringLong2DigitTrie regionsPerTaxid) {
        return new MyFastaReader(intConfigValue(GSConfigKey.FASTA_LINE_SIZE_BYTES),
                relevantNodes,
                isIncludeRefSeqFna() ? accessionMapGoal.get() : null,
                intConfigValue(GSConfigKey.KMER_SIZE),
                intConfigValue(GSConfigKey.MAX_GENOMES_PER_TAXID),
                (Rank) configValue(GSConfigKey.MAX_GENOMES_PER_TAXID_RANK),
                longConfigValue(GSConfigKey.MAX_KMERS_PER_TAXID),
                intConfigValue(GSConfigKey.MAX_DUST),
                intConfigValue(GSConfigKey.STEP_SIZE),
                booleanConfigValue(GSConfigKey.COMPLETE_GENOMES_ONLY),
                regionsPerTaxid,
                booleanConfigValue(GSConfigKey.ENABLE_LOWERCASE_BASES));
    }

    protected class MyFastaReader extends AbstractStoreFastaReader {
        private SmallTaxTree.SmallTaxIdNode smallNode;

        public MyFastaReader(int bufferSize, Set<TaxTree.TaxIdNode> taxNodes, AccessionMap accessionMap,
                             int k, int maxGenomesPerTaxId, Rank maxGenomesPerTaxIdRank, long maxKmersPerTaxId, int maxDust, int stepSize, boolean completeGenomesOnly, StringLong2DigitTrie regionsPerTaxid, boolean enableLowerCaseBases) {
            super(bufferSize, taxNodes, accessionMap, k, maxGenomesPerTaxId, maxGenomesPerTaxIdRank, maxKmersPerTaxId, maxDust, stepSize, completeGenomesOnly, regionsPerTaxid, enableLowerCaseBases);
        }

        // The following 3 override ensure that all k-mers for considers tax nodes are used
        // so that the clustering phase can work with the maximum amount of data.
        @Override
        protected void infoLine() {
            if (ignoreMap) {
                node = mappedNode;
            } else {
                updateNodeFromInfoLine();
            }

            if (node != null && taxNodes.contains(node)) {
                includeRegion = true;
                // This can become null if the node is not in the database:
                smallNode = smallTaxTree.getNodeByTaxId(node.getTaxId());
            }
        }

        @Override
        protected void endRegion() {
            // Intentionally empty.
        }

        @Override
        public boolean isAllowMoreKmers() {
            return true;
        }

        @Override
        protected boolean handleStore() {
            long kmer = byteRingBuffer.getStandardKMer();
            SmallTaxTree.SmallTaxIdNode storedNode = kMerSortedArray.getLong(kmer, null);
            if (storedNode != null) {
                if (FTConfigKey.RefinementPosition.getMatchingNodeFor(storedNode, refinementIntervals) != null) {
                    int index = smallNode == null ? OTHER_VALUE : smallNode.storeIndex;
                    if (!filter.containsLongInt(kmer, index)) {
                        if (multiThreading) {
                            synchronized (filter) {
                                // This is a trick to enable more parallelism -
                                // check again after synchronized to avoid synchronized further outside...
                                if (!filter.containsLongInt(kmer, index)) {
                                    filter.putLongInt(kmer, index);
                                    return true;
                                }
                            }
                        } else {
                            filter.putLongInt(kmer, index);
                            return true;
                        }
                    }
                }
            }
            return false;
        }
    }
}
