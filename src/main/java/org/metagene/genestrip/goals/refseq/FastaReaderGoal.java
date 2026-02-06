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
package org.metagene.genestrip.goals.refseq;

import me.tongfei.progressbar.ProgressBar;
import org.metagene.genestrip.ExecutionContext;
import org.metagene.genestrip.GSConfigKey;
import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.fastq.AbstractFastqReader;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.GoalKey;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.refseq.AbstractRefSeqFastaReader;
import org.metagene.genestrip.refseq.AbstractStoreFastaReader;
import org.metagene.genestrip.refseq.RefSeqCategory;
import org.metagene.genestrip.tax.TaxTree;
import org.metagene.genestrip.util.SimpleBlockinqQueue;
import org.metagene.genestrip.util.StringLongDigitTrie;
import org.metagene.genestrip.util.progressbar.GSProgressBarCreator;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;

public abstract class FastaReaderGoal<T, P extends GSProject> extends ObjectGoal<T, P> {
    protected final ObjectGoal<Set<RefSeqCategory>, P> categoriesGoal;
    protected final ObjectGoal<Set<TaxTree.TaxIdNode>, P> taxNodesGoal;
    protected final RefSeqFnaFilesDownloadGoal fnaFilesGoal;
    protected final ObjectGoal<Map<File, TaxTree.TaxIdNode>, P> additionalGoal;

    private final ExecutionContext bundle;

    private boolean dump;
    private int doneCounter;
    private ProgressBar progressBar;

    public FastaReaderGoal(P project, GoalKey key, ExecutionContext bundle, ObjectGoal<Set<RefSeqCategory>, P> categoriesGoal,
                           ObjectGoal<Set<TaxTree.TaxIdNode>, P> taxNodesGoal, RefSeqFnaFilesDownloadGoal fnaFilesGoal,
                           ObjectGoal<Map<File, TaxTree.TaxIdNode>, P> additionalGoal, Goal<P>... dependencies) {
        super(project, key, Goal.append(dependencies, categoriesGoal, taxNodesGoal, fnaFilesGoal, additionalGoal));
        this.categoriesGoal = categoriesGoal;
        this.taxNodesGoal = taxNodesGoal;
        this.fnaFilesGoal = fnaFilesGoal;
        this.additionalGoal = additionalGoal;
        this.bundle = bundle;
    }

    public void readFastas() throws IOException {
        BlockingQueue<FileAndNode> blockingQueue = null;
        AbstractRefSeqFastaReader.StringLong2DigitTrie regionsPerTaxid = new AbstractRefSeqFastaReader.StringLong2DigitTrie();
        // For minUpdate the fna files must be read in the same order as during the goals from before.
        // Otherwise, the wrong k-mers will be compared to the ones from the DB.
        // For !minUpdate multi-threading can be enabled.
        if (bundle.getThreads() > 0) {
            blockingQueue = createBlockingQueue(intConfigValue(GSConfigKey.THREAD_QUEUE_SIZE));
            for (int i = 0; i < bundle.getThreads(); i++) {
                bundle.execute(createFastaReaderRunnable(i, blockingQueue, regionsPerTaxid));
            }
        }
        AbstractRefSeqFastaReader fastaReader = createFastaReader(regionsPerTaxid);

        int sumFiles = 0;
        List<File> refSeqFiles = isIncludeRefSeqFna() ? fnaFilesGoal.getFiles() : Collections.emptyList();
        sumFiles += refSeqFiles.size();
        Map<File, TaxTree.TaxIdNode> additionalMap = additionalGoal == null ? null : additionalGoal.get();
        sumFiles += additionalMap == null ? 0 : additionalMap.size();
        try (ProgressBar pb = (progressBar = createProgressBar(sumFiles))) {
            doneCounter = 0;
            for (File fnaFile : refSeqFiles) {
                RefSeqCategory cat = fnaFilesGoal.getCategoryForFile(fnaFile);
                if (categoriesGoal.get().contains(cat)) {
                    if (blockingQueue == null) {
                        fastaReader.readFasta(fnaFile);
                    } else {
                        try {
                            doneCounter++;
                            blockingQueue.put(new DBGoal.FileAndNode(fnaFile, null));
                        } catch (InterruptedException e) {
                            throw new RuntimeException(e);
                        }
                    }
                }
                checkAndLogConsumerThreadProblem();
            }
            if (additionalMap != null) {
                for (File additionalFasta : additionalMap.keySet()) {
                    if (blockingQueue == null) {
                        fastaReader.ignoreAccessionMap(additionalMap.get(additionalFasta));
                        fastaReader.readFasta(additionalFasta);
                    } else {
                        try {
                            doneCounter++;
                            blockingQueue.put(new DBGoal.FileAndNode(additionalFasta, additionalMap.get(additionalFasta)));
                        } catch (InterruptedException e) {
                            throw new RuntimeException(e);
                        }
                    }
                    checkAndLogConsumerThreadProblem();
                }
            }
            // Gentle polling and waiting until all consumers are done.
            while (doneCounter > 0 && !dump) {
                checkAndLogConsumerThreadProblem();
                try {
                    Thread.sleep(100);
                } catch (InterruptedException e) {
                    // Ignore.
                }
            }
        }
        bundle.clearThrowableList();
        afterReadFastas(regionsPerTaxid);
    }

    protected BlockingQueue<FileAndNode> createBlockingQueue(int maxQueueSize) {
        // This simple blocking queue performs better than ArrayBlockingQueue.
        return new SimpleBlockinqQueue<>(maxQueueSize);
        //return new ArrayBlockingQueue<>(maxQueueSize);
    }

    protected boolean isIncludeRefSeqFna() {
        return booleanConfigValue(GSConfigKey.REF_SEQ_DB);
    }

    protected void afterReadFastas(AbstractRefSeqFastaReader.StringLong2DigitTrie regionsPerTaxid) {
    }

    protected ProgressBar createProgressBar(int max) {
        return booleanConfigValue(GSConfigKey.PROGRESS_BAR) ?
                GSProgressBarCreator.newGSProgressBar(getKey().getName(), max, 60000, " files", null, getLogger(), false) :
                null;
    }

    protected void checkAndLogConsumerThreadProblem() {
        if (!bundle.getThrowableList().isEmpty()) {
            for (Throwable t : bundle.getThrowableList()) {
                if (getLogger().isErrorEnabled()) {
                    getLogger().error("Error in consumer thread: ", t);
                }
            }
            bundle.clearThrowableList();
            throw new RuntimeException("Error(s) in consumer thread(s).");
        }
    }

    protected Runnable createFastaReaderRunnable(int i,
                                                 BlockingQueue<DBGoal.FileAndNode> blockingQueue,
                                                 AbstractRefSeqFastaReader.StringLong2DigitTrie regionsPerTaxid) {
        AbstractRefSeqFastaReader fastaReader = createFastaReader(regionsPerTaxid);
        return new Runnable() {
            @Override
            public void run() {
                while (!dump) {
                    try {
                        try {
                            DBGoal.FileAndNode fileAndNode = blockingQueue.take();
                            fastaReader.ignoreAccessionMap(fileAndNode.getNode());
                            fastaReader.readFasta(fileAndNode.getFile());
                            if (progressBar != null) {
                                progressBar.step();
                            }
                        } finally {
                            doneCounter--;
                        }
                    } catch (IOException e) {
                        throw new RuntimeException(e);
                    } catch (InterruptedException e) {
                        if (!dump) {
                            throw new RuntimeException(e);
                        }
                    }
                }
            }
        };
    }

    protected abstract AbstractRefSeqFastaReader createFastaReader(AbstractRefSeqFastaReader.StringLong2DigitTrie regionsPerTaxid);

    public void dump() {
        super.dump();
        cleanUpThreads();
    }

    protected void cleanUpThreads() {
        dump = true;
        bundle.interruptAll();
    }

    protected static final class FileAndNode {
        private final File file;
        private final TaxTree.TaxIdNode node;

        public FileAndNode(File file, TaxTree.TaxIdNode node) {
            this.file = file;
            this.node = node;
        }

        public File getFile() {
            return file;
        }

        public TaxTree.TaxIdNode getNode() {
            return node;
        }
    }
}
