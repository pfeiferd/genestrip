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
import org.metagene.genestrip.GSProject.FileType;
import org.metagene.genestrip.goals.MultiFileGoal;
import org.metagene.genestrip.io.StreamProvider;
import org.metagene.genestrip.io.StreamingFileResource;
import org.metagene.genestrip.io.StreamingResource;
import org.metagene.genestrip.io.StreamingResourceStream;
import org.metagene.genestrip.kraken.KrakenExecutor;
import org.metagene.genestrip.kraken.KrakenResultListener;
import org.metagene.genestrip.kraken.KrakenResultProcessor;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.match.MatchingResult;
import org.metagene.genestrip.store.Database;
import org.metagene.genestrip.store.KMerUniqueCounter;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;
import org.metagene.genestrip.util.DigitTrie;

public class KrakenResCountGoal extends ObjectGoal<Map<String, List<KrakenResCountGoal.KrakenResStats>>, GSProject> {
    protected final ObjectGoal<Map<String, StreamingResourceStream>, GSProject> fastqMapGoal;
    private final ObjectGoal<Set<TaxIdNode>, GSProject> taxNodesGoal;

    @SafeVarargs
    public KrakenResCountGoal(GSProject project,
                              ObjectGoal<Map<String, StreamingResourceStream>, GSProject> fastqMapGoal,
                              ObjectGoal<Set<TaxIdNode>, GSProject> taxNodesGoal, Goal<GSProject>... deps) {
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
                File filteredFile = null;
                File krakenOutStyleFile = null;
                StreamingResourceStream fastqs = map.get(key);
                List<File> files = new ArrayList<>();
                for (StreamingResource s : fastqs) {
                    files.add(((StreamingFileResource) s).getFile());
                }
                countResults.put(key, computeStats(key, files));
            }
            set(countResults);
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

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
                                afterReadMatch(krakenTaxid, readDescriptor);
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

            if (getLogger().isInfoEnabled()) {
                getLogger().info("Finished kraken");
            }

            List<KrakenResStats> list = new ArrayList<KrakenResStats>();
            countingTrie.collect(list);
            return list;
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    protected void afterReadMatch(String krakenTaxid, byte[] readDescriptor) {
    }

    public static class KrakenResStats implements Comparable<KrakenResStats>, Serializable {
        private static final long serialVersionUID = 1L;

        private String taxid;

        private long reads;
        private long kmers;
        private long kmersInMatchingReads;

        @Override
        public int compareTo(KrakenResStats o) {
            return taxid.compareTo(o.taxid);
        }

        public long getReads() {
            return reads;
        }

        public long getKmers() {
            return kmers;
        }

        public long getKmersInMatchingReads() {
            return kmersInMatchingReads;
        }

        public String getTaxid() {
            return taxid;
        }
    }
}
