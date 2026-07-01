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
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.GoalKey;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.refseq.AbstractRefSeqFastaReader;
import org.metagene.genestrip.refseq.RefSeqCategory;
import org.metagene.genestrip.tax.TaxTree;
import org.metagene.genestrip.util.SimpleBlockingQueue;
import org.metagene.genestrip.util.progressbar.GSProgressBarCreator;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.atomic.AtomicInteger;

/**
 * Abstract base for goals that read the downloaded RefSeq {@code .fna} files (and any additional
 * FASTA files) region by region, optionally in parallel via a pool of consumer threads, dispatching
 * each region to a subclass-provided {@link AbstractRefSeqFastaReader}.
 *
 * @param <T> the type of result produced by this goal
 * @param <P> the project type
 */
public abstract class FastaReaderGoal<T, P extends GSProject> extends ObjectGoal<T, P> {
    /** The goal supplying the set of RefSeq categories to read. */
    protected final ObjectGoal<Set<RefSeqCategory>, P> categoriesGoal;
    /** The goal supplying the set of required taxonomy nodes. */
    protected final ObjectGoal<Set<TaxTree.TaxIdNode>, P> taxNodesGoal;
    /** The goal supplying the downloaded RefSeq {@code .fna} files. */
    protected final RefSeqFnaFilesDownloadGoal fnaFilesGoal;
    /** The goal supplying additional FASTA files mapped to their tax node. */
    protected final ObjectGoal<Map<File, TaxTree.TaxIdNode>, P> additionalGoal;

    private final ExecutionContext bundle;

    // volatile / atomic: written by dump() and the consumer threads, read by the producer's spin loop.
    private volatile boolean dump;
    private final AtomicInteger doneCounter = new AtomicInteger();
    private ProgressBar progressBar;

    /**
     * Creates the reader goal with its category, tax-node, RefSeq-file-download and additional-file
     * dependency goals and the execution context that supplies the worker threads.
     *
     * @param project        the project this goal belongs to
     * @param key            the key identifying this goal
     * @param bundle         the execution context supplying the worker threads
     * @param categoriesGoal the goal supplying the RefSeq categories to read
     * @param taxNodesGoal   the goal supplying the required taxonomy nodes
     * @param fnaFilesGoal   the goal supplying the downloaded RefSeq {@code .fna} files
     * @param additionalGoal the goal supplying additional FASTA files mapped to their tax node
     * @param dependencies   any further goals this goal depends on
     */
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

    /**
     * Reads all relevant RefSeq FASTA files and any additional FASTA files, single-threaded or via
     * the configured pool of consumer threads, then invokes {@link #afterReadFastas}.
     *
     * @throws IOException if reading a FASTA file fails
     */
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
            doneCounter.set(0);
            for (File fnaFile : refSeqFiles) {
                RefSeqCategory cat = fnaFilesGoal.getCategoryForFile(fnaFile);
                if (categoriesGoal.get().contains(cat)) {
                    if (blockingQueue == null) {
                        fastaReader.readFasta(fnaFile);
                    } else {
                        try {
                            doneCounter.incrementAndGet();
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
                            doneCounter.incrementAndGet();
                            blockingQueue.put(new DBGoal.FileAndNode(additionalFasta, additionalMap.get(additionalFasta)));
                        } catch (InterruptedException e) {
                            throw new RuntimeException(e);
                        }
                    }
                    checkAndLogConsumerThreadProblem();
                }
            }
            // Gentle polling and waiting until all consumers are done.
            while (doneCounter.get() > 0 && !dump) {
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

    /**
     * Creates the bounded queue that feeds FASTA files to the consumer threads.
     *
     * @param maxQueueSize the maximum number of queued files
     * @return the newly created blocking queue
     */
    protected BlockingQueue<FileAndNode> createBlockingQueue(int maxQueueSize) {
        // This simple blocking queue performs better than ArrayBlockingQueue.
        return new SimpleBlockingQueue<>(maxQueueSize);
        //return new ArrayBlockingQueue<>(maxQueueSize);
    }

    /**
     * Whether the downloaded RefSeq {@code .fna} files are part of this goal's input.
     *
     * @return {@code true} if the RefSeq {@code .fna} files should be read
     */
    protected boolean isIncludeRefSeqFna() {
        return booleanConfigValue(GSConfigKey.REF_SEQ_DB);
    }

    /**
     * Hook invoked once all FASTA files have been read; the default implementation does nothing.
     *
     * @param regionsPerTaxid the trie tracking how many regions were seen per taxid
     */
    protected void afterReadFastas(AbstractRefSeqFastaReader.StringLong2DigitTrie regionsPerTaxid) {
    }

    /**
     * Creates the progress bar spanning the given number of files, or {@code null} if progress bars
     * are disabled.
     *
     * @param max the number of files the progress bar spans
     * @return the progress bar, or {@code null} if progress bars are disabled
     */
    protected ProgressBar createProgressBar(int max) {
        return booleanConfigValue(GSConfigKey.PROGRESS_BAR) ?
                GSProgressBarCreator.newGSProgressBar(getKey().getName(), max, 60000, " files", null, getLogger(), false) :
                null;
    }

    /**
     * Logs and rethrows any exceptions collected from the consumer threads.
     */
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

    /**
     * Creates a consumer {@link Runnable} that takes files from the queue and reads them with its
     * own FASTA reader.
     *
     * @param i               the index of the consumer thread
     * @param blockingQueue   the queue supplying files to read
     * @param regionsPerTaxid the shared trie tracking how many regions were seen per taxid
     * @return the consumer runnable
     */
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
                            doneCounter.decrementAndGet();
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

    /**
     * Creates the FASTA reader that processes each region; called once per reader thread. The shared
     * {@code regionsPerTaxid} trie tracks how many regions have been seen per taxid.
     *
     * @param regionsPerTaxid the shared trie tracking how many regions were seen per taxid
     * @return the FASTA reader that processes each region
     */
    protected abstract AbstractRefSeqFastaReader createFastaReader(AbstractRefSeqFastaReader.StringLong2DigitTrie regionsPerTaxid);

    /**
     * In addition to discarding the result, stops and interrupts the consumer threads.
     */
    public void dump() {
        super.dump();
        cleanUpThreads();
    }

    /**
     * Signals the consumer threads to stop and interrupts any that are blocked.
     */
    protected void cleanUpThreads() {
        dump = true;
        bundle.interruptAll();
    }

    /**
     * A FASTA file paired with an optional tax node; a non-null node marks an additional FASTA whose
     * node bypasses the accession map.
     */
    protected static final class FileAndNode {
        private final File file;
        private final TaxTree.TaxIdNode node;

        /**
         * Creates a pairing of the given FASTA file and optional tax node.
         *
         * @param file the FASTA file
         * @param node the associated tax node, or {@code null}
         */
        public FileAndNode(File file, TaxTree.TaxIdNode node) {
            this.file = file;
            this.node = node;
        }

        /**
         * Returns the FASTA file.
         *
         * @return the file
         */
        public File getFile() {
            return file;
        }

        /**
         * Returns the associated tax node, or {@code null} if there is none.
         *
         * @return the tax node
         */
        public TaxTree.TaxIdNode getNode() {
            return node;
        }
    }
}
