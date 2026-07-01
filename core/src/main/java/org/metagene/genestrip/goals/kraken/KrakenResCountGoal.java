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
package org.metagene.genestrip.goals.kraken;

import java.io.*;
import java.util.*;

import org.metagene.genestrip.GSConfigKey;
import org.metagene.genestrip.GSGoalKey;
import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.io.StreamingFileResource;
import org.metagene.genestrip.io.StreamingResource;
import org.metagene.genestrip.io.StreamingResourceStream;
import org.metagene.genestrip.kraken.KrakenExecutor;
import org.metagene.genestrip.kraken.KrakenResultListener;
import org.metagene.genestrip.kraken.KrakenResultProcessor;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;
import org.metagene.genestrip.util.DigitTrie;

/**
 * Runs Kraken on each input fastq set and produces per-key, per-taxid read and k-mer statistics
 * ({@link KrakenResStats}), optionally restricted to a given set of tax ids.
 *
 * @param <P> the project type
 */
public class KrakenResCountGoal<P extends GSProject> extends ObjectGoal<Map<String, List<KrakenResCountGoal.KrakenResStats>>, P> {
    /** Goal supplying the per-key map of fastq input resources to run Kraken on. */
    protected final ObjectGoal<Map<String, StreamingResourceStream>, P> fastqMapGoal;
    private final ObjectGoal<Set<TaxIdNode>, P> taxNodesGoal;
    private AfterMatchCallback afterMatchCallback;

    /**
     * Creates the goal.
     *
     * @param project      the project this goal belongs to
     * @param fastqMapGoal goal supplying the per-key map of fastq input resources
     * @param taxNodesGoal goal supplying the tax id nodes to restrict counting to, may be {@code null}
     * @param deps         additional goals this goal depends on
     */
    @SafeVarargs
    public KrakenResCountGoal(P project,
                              ObjectGoal<Map<String, StreamingResourceStream>, P> fastqMapGoal,
                              ObjectGoal<Set<TaxIdNode>, P> taxNodesGoal, Goal<P>... deps) {
        super(project, GSGoalKey.KRAKENCOUNT, append(deps, taxNodesGoal, fastqMapGoal));
        this.fastqMapGoal = fastqMapGoal;
        this.taxNodesGoal = taxNodesGoal;
    }

    @Override
    protected void doMakeThis() {
        try {
            Map<String, List<KrakenResStats>> countResults = new HashMap<>();
            Map<String, StreamingResourceStream> map = fastqMapGoal.get();

            for (String key : map.keySet()) {
                StreamingResourceStream fastqs = map.get(key);
                List<File> files = new ArrayList<>();
                for (StreamingResource s : fastqs) {
                    files.add(((StreamingFileResource) s).getFile());
                }
                List<KrakenResStats> res = computeStats(key, files);
                countResults.put(key, res);
                if (afterMatchCallback != null) {
                    afterMatchCallback.afterKey(key, res);
                }
            }
            set(countResults);
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    /**
     * Runs Kraken on the given fastq files and returns the per-taxid statistics, restricted to the configured
     * tax ids when a tax-nodes goal is set.
     *
     * @param key    the input key the statistics are computed for
     * @param fastqs the fastq files to run Kraken on
     * @return the per-taxid Kraken statistics
     * @throws IOException if running Kraken or reading its output fails
     */
    protected List<KrakenResStats> computeStats(String key, List<File> fastqs) throws IOException {
        DigitTrie<KrakenResStats> countingTrie = new DigitTrie<KrakenResStats>() {
            @Override
            protected KrakenResStats createInGet(String taxid, Object createContext) {
                KrakenResStats stats = new KrakenResStats();
                stats.taxid = taxid;
                return stats;
            }
        };

        final Set<String> taxIds;

        if (taxNodesGoal != null) {
            taxIds = new HashSet<String>();
            for (TaxIdNode node : taxNodesGoal.get()) {
                taxIds.add(node.getTaxId());
            }
        } else {
            taxIds = null;
        }

        KrakenExecutor krakenExecutor = new KrakenExecutor(stringConfigValue(GSConfigKey.KRAKEN_BIN),
                stringConfigValue(GSConfigKey.KRAKEN_EXEC_EXPR)) {
            @Override
            protected void handleOutputStream(InputStream stream, OutputStream out) throws IOException {
                KrakenResultProcessor parser = new KrakenResultProcessor(65536);

                parser.process(new BufferedInputStream(stream), new KrakenResultListener() {
                    private long lastLine = -1;

                    @Override
                    public void newTaxIdForRead(long lineCount, byte[] readDescriptor, String krakenTaxid, int bps,
                                                int pos, String kmerTaxid, int hitLength, byte[] output) {
                        if (taxIds == null || taxIds.contains(kmerTaxid)) {
                            KrakenResStats stats = countingTrie.get(kmerTaxid, true);
                            stats.kmers += hitLength;
                        }
                        if (lineCount != lastLine) {
                            lastLine = lineCount;
                            if (taxIds == null || taxIds.contains(krakenTaxid)) {
                                KrakenResStats stats = countingTrie.get(krakenTaxid, true);
                                stats.reads++;
                                if (kmerTaxid != null && kmerTaxid.equals(krakenTaxid)) {
                                    stats.kmersInMatchingReads += hitLength;
                                }
                                if (afterMatchCallback != null) {
                                    afterMatchCallback.afterMatch(krakenTaxid, readDescriptor);
                                }
                            }
                        }
                    }
                });
            }
        };
        if (getLogger().isInfoEnabled()) {
            String execLine = krakenExecutor.genExecLine(stringConfigValue(GSConfigKey.KRAKEN_DB), fastqs, null);
            getLogger().info("Run kraken with " + execLine);
        }
        try {
            if (krakenExecutor.isWithFileForOutput()) {
                throw new IOException("This goal does not work with an outfile as a parameter (like in krakenuniq)");
            }
            /*
            File outFile = getProject().getOutputFile(GSGoalKey.KRAKENCOUNT.getName(), key, null, FileType.KRAKEN_OUT_RES, true);
            if (getLogger().isInfoEnabled()) {
                getLogger().info("Writing kraken out to: " + outFile.getAbsolutePath());
            }
            try (OutputStream out = new FileOutputStream(outFile)) {
             */
                krakenExecutor.execute2(stringConfigValue(GSConfigKey.KRAKEN_DB), fastqs, null, null, System.err);
            //}

            if (getLogger().isDebugEnabled()) {
                getLogger().debug("Finished kraken");
            }

            List<KrakenResStats> list = new ArrayList<KrakenResStats>();
            countingTrie.collect(list);
            return list;
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    /**
     * Hook invoked for each classified read; does nothing by default.
     *
     * @param krakenTaxid    the tax id Kraken assigned to the read
     * @param readDescriptor the read descriptor bytes
     */
    protected void afterReadMatch(String krakenTaxid, byte[] readDescriptor) {
    }

    /**
     * Per-taxid Kraken statistics: number of reads, matched k-mers, and k-mers within matching reads.
     */
    public static class KrakenResStats implements Comparable<KrakenResStats>, Serializable {
        private static final long serialVersionUID = 1L;

        /** The tax id these statistics belong to. */
        private String taxid;

        /** The number of reads classified to the tax id. */
        private long reads;
        /** The number of matched k-mers counted for the tax id. */
        private long kmers;
        /** The number of k-mers within reads matching the tax id. */
        private long kmersInMatchingReads;

        /**
         * Creates an empty statistics record.
         */
        public KrakenResStats() {
        }

        @Override
        public int compareTo(KrakenResStats o) {
            return taxid.compareTo(o.taxid);
        }

        /**
         * Returns the number of reads classified to the tax id.
         *
         * @return the number of reads classified to the tax id
         */
        public long getReads() {
            return reads;
        }

        /**
         * Returns the number of matched k-mers counted for the tax id.
         *
         * @return the number of matched k-mers counted for the tax id
         */
        public long getKmers() {
            return kmers;
        }

        /**
         * Returns the number of k-mers within reads matching the tax id.
         *
         * @return the number of k-mers within reads matching the tax id
         */
        public long getKmersInMatchingReads() {
            return kmersInMatchingReads;
        }

        /**
         * Returns the tax id these statistics belong to.
         *
         * @return the tax id these statistics belong to
         */
        public String getTaxid() {
            return taxid;
        }
    }

    /**
     * Sets the callback invoked while counting Kraken results.
     *
     * @param afterMatchCallback the callback to invoke, may be {@code null}
     */
    public void setAfterMatchCallback(AfterMatchCallback afterMatchCallback) {
        this.afterMatchCallback = afterMatchCallback;
    }

    /**
     * Callback invoked while counting Kraken results: once per finished key and once per classified read.
     */
    public interface AfterMatchCallback {
        /**
         * Invoked once the statistics for a key have been computed.
         *
         * @param key the input key that was finished
         * @param res the per-taxid statistics computed for the key
         */
        public void afterKey(String key, List<KrakenResStats> res);

        /**
         * Invoked for each classified read.
         *
         * @param krakenTaxid    the tax id Kraken assigned to the read
         * @param readDescriptor the read descriptor bytes
         */
        public void afterMatch(String krakenTaxid, byte[] readDescriptor);
    }
}
