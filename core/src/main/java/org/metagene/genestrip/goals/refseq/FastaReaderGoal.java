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
import org.metagene.genestrip.tax.TaxNodeSelection;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.atomic.AtomicInteger;

/**
 * Abstract base for goals that read the downloaded RefSeq {@code .fna} files (and any additional
 * FASTA files) contig by contig, optionally in parallel via a pool of consumer threads, dispatching
 * each contig to a subclass-provided {@link AbstractRefSeqFastaReader}.
 *
 * @param <T> the type of result produced by this goal
 * @param <P> the project type
 */
public abstract class FastaReaderGoal<T, P extends GSProject> extends ObjectGoal<T, P> {
    /** The goal supplying the set of RefSeq categories to read. */
    protected final ObjectGoal<Set<RefSeqCategory>, P> categoriesGoal;
    /** The goal supplying the set of required taxonomy nodes. */
    protected final ObjectGoal<TaxNodeSelection, P> taxNodesGoal;
    /** The goal supplying the downloaded RefSeq {@code .fna} files. */
    protected final RefSeqFnaFilesDownloadGoal fnaFilesGoal;
    /** The goal supplying additional FASTA files mapped to their tax node. */
    protected final ObjectGoal<Map<File, TaxTree.TaxIdNode>, P> additionalGoal;

    private final ExecutionContext bundle;

    /**
     * Queued once per consumer when a pass has read everything it was given, telling each to return;
     * see {@link #endOfPass(BlockingQueue)}. Recognised by identity, so it carries neither file nor
     * node.
     */
    private static final FileAndNode END_OF_PASS = new FileAndNode(null, null);

    // volatile / atomic: written by dump() and the consumer threads, read by the producer's spin loop.
    private volatile boolean dump;
    private final AtomicInteger doneCounter = new AtomicInteger();
    private ProgressBar progressBar;
    // Decided once, when the goal is created, because it also decides whether the RefSeq download
    // goal becomes a dependency - and dependencies are fixed at construction time.
    private final boolean includeRefSeqFna;

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
                           ObjectGoal<TaxNodeSelection, P> taxNodesGoal, RefSeqFnaFilesDownloadGoal fnaFilesGoal,
                           ObjectGoal<Map<File, TaxTree.TaxIdNode>, P> additionalGoal, Goal<P>... dependencies) {
        this(project, key, bundle, categoriesGoal, taxNodesGoal, fnaFilesGoal, additionalGoal,
                project.booleanConfigValue(GSConfigKey.REF_SEQ_DB), dependencies);
    }

    /**
     * Creates the reader goal, stating explicitly whether it reads the RefSeq release.
     * <p>
     * A goal that does not read it does not depend on it either: the download goal is then left out
     * of the dependencies, so that requesting such a goal does not fetch and verify a RefSeq release
     * whose content it will never look at. Only the goal that computes the lowest common ancestors
     * passes {@code true} regardless of the configuration - see {@code DBGoal} for why.
     *
     * @param project           the project this goal belongs to
     * @param key               the key identifying this goal
     * @param bundle            the execution context supplying the worker threads
     * @param categoriesGoal    the goal supplying the RefSeq categories to read
     * @param taxNodesGoal      the goal supplying the required taxonomy nodes
     * @param fnaFilesGoal      the goal supplying the downloaded RefSeq {@code .fna} files
     * @param additionalGoal    the goal supplying additional FASTA files mapped to their tax node
     * @param includeRefSeqFna  whether this goal reads the RefSeq release at all
     * @param dependencies      any further goals this goal depends on
     */
    public FastaReaderGoal(P project, GoalKey key, ExecutionContext bundle, ObjectGoal<Set<RefSeqCategory>, P> categoriesGoal,
                           ObjectGoal<TaxNodeSelection, P> taxNodesGoal, RefSeqFnaFilesDownloadGoal fnaFilesGoal,
                           ObjectGoal<Map<File, TaxTree.TaxIdNode>, P> additionalGoal, boolean includeRefSeqFna,
                           Goal<P>... dependencies) {
        super(project, key, Goal.append(dependencies, categoriesGoal, taxNodesGoal,
                includeRefSeqFna ? fnaFilesGoal : null, additionalGoal));
        this.categoriesGoal = categoriesGoal;
        this.taxNodesGoal = taxNodesGoal;
        this.fnaFilesGoal = fnaFilesGoal;
        this.additionalGoal = additionalGoal;
        this.bundle = bundle;
        this.includeRefSeqFna = includeRefSeqFna;
    }

    /**
     * Reads all relevant RefSeq FASTA files and any additional FASTA files, single-threaded or via
     * the configured pool of consumer threads, then invokes {@link #afterReadFastas}.
     *
     * @throws IOException if reading a FASTA file fails
     */
    public void readFastas() throws IOException {
        // Cleared here so that a goal aborted through dump() is not finished for good; a pass that
        // ends normally leaves it false anyway, because it ends its own consumers instead - see
        // endOfPass(). It did not always: the flag doubled as the end-of-pass signal, every subclass
        // set it from its own finally, and nothing put it back, so a goal that had read once queued
        // everything the next time and read none of it while reporting success.
        readyForAnotherPass();
        BlockingQueue<FileAndNode> blockingQueue = null;
        AbstractRefSeqFastaReader.StringLong2DigitTrie contigsPerTaxid = new AbstractRefSeqFastaReader.StringLong2DigitTrie();
        // Reading order does not affect what a k-mer ends up mapped to, whatever `updateScope' is: the
        // update only touches k-mers already present in the DB, and each is merged into the lowest
        // common ancestor of its stored node and the nodes of the contigs it occurs in. That merge is
        // commutative and associative, so the outcome depends only on the set of included contigs and
        // not on the order in which the fna files are read.
        //
        // Which contigs those are is a different matter, and is not order-independent once
        // maxGenomesPerTaxid binds: a genome is admitted while its taxon's count is still below the
        // limit, and which thread reaches a taxon first decides which of its genomes get in. Two
        // builds of the same project can therefore hold different genomes. Within one build they
        // cannot: the sizing pass makes the selection, freezes it, and every pass after it follows
        // that one set - which is what keeps the fill from filing k-mers at a node the counting pass
        // never registered. At the default, where the limit does not bind, the set of included
        // contigs is fixed and everything here is exact.
        if (bundle.getThreads() > 0) {
            blockingQueue = createBlockingQueue(intConfigValue(GSConfigKey.THREAD_QUEUE_SIZE));
            for (int i = 0; i < bundle.getThreads(); i++) {
                bundle.execute(createFastaReaderRunnable(i, blockingQueue, contigsPerTaxid));
            }
        }
        AbstractRefSeqFastaReader fastaReader = createFastaReader(contigsPerTaxid);

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
                        // Stated rather than left to the reader's initial state: a release file is
                        // the one case where no node comes with the file, and a pass that tells the
                        // two apart - see AbstractRefSeqFastaReader#isRefSeqReleaseContig() - must
                        // not depend on the release happening to be read before the project's own
                        // fastas. The queued path below passes the same null through FileAndNode.
                        fastaReader.ignoreAccessionMap(null);
                        fastaReader.readFasta(fnaFile);
                    } else {
                        try {
                            doneCounter.incrementAndGet();
                            blockingQueue.put(new FileAndNode(fnaFile, null));
                        } catch (InterruptedException e) {
                            // A consumer that dies records its throwable and calls interruptAll(),
                            // which interrupts this thread too - so an interrupt here is usually the
                            // symptom and never the cause. Reporting the cause first keeps the real
                            // stack trace from being replaced by an InterruptedException nobody can
                            // act on. If no consumer failed, the interrupt stands on its own.
                            checkAndLogConsumerThreadProblem();
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
                            blockingQueue.put(new FileAndNode(additionalFasta, additionalMap.get(additionalFasta)));
                        } catch (InterruptedException e) {
                            // A consumer that dies records its throwable and calls interruptAll(),
                            // which interrupts this thread too - so an interrupt here is usually the
                            // symptom and never the cause. Reporting the cause first keeps the real
                            // stack trace from being replaced by an InterruptedException nobody can
                            // act on. If no consumer failed, the interrupt stands on its own.
                            checkAndLogConsumerThreadProblem();
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
            endOfPass(blockingQueue);
        }
        bundle.clearThrowableList();
        afterReadFastas(contigsPerTaxid);
    }

    /**
     * Ends the consumer threads of a pass that has read everything it was given.
     * <p>
     * They do not end on their own: a consumer loops on {@code blockingQueue.take()} and the only
     * thing that ever released it used to be the flag {@link #dump}, which is meant for aborting.
     * Leaving a pass unended is not merely untidy - the pool is a fixed one, sized to exactly as many
     * threads as a pass has consumers, so every finished pass would keep all of them parked and the
     * next pass would get no thread at all. So each subclass ended the pass itself, by setting the
     * abort flag from its own {@code finally}, and since nothing put that flag back the goal could
     * never read a second time. Ending the pass here rather than there is what lets the flag mean
     * only what its name says.
     * <p>
     * One sentinel per consumer, rather than the flag: it is put only after the wait loop above, so
     * no file is ever left queued behind it, and each consumer takes exactly one and returns. No
     * interrupt is involved and no window exists in which the flag has to hold a particular value,
     * which is what a flag shared between two passes could not offer.
     *
     * @param blockingQueue the queue the consumers of this pass are waiting on, or {@code null} when
     *                      the pass read single-threaded and has no consumers
     */
    private void endOfPass(BlockingQueue<FileAndNode> blockingQueue) {
        if (blockingQueue == null || dump) {
            // An aborted pass has ended its consumers through the flag already.
            return;
        }
        try {
            for (int i = 0; i < bundle.getThreads(); i++) {
                blockingQueue.put(END_OF_PASS);
            }
        } catch (InterruptedException e) {
            // The reading is complete and its result is in hand, so this is not a failure of the pass.
            // Being interrupted here means someone is aborting, and the flag they set ends the
            // consumers anyway.
            Thread.currentThread().interrupt();
        }
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
    protected final boolean isIncludeRefSeqFna() {
        return includeRefSeqFna;
    }

    /**
     * Hook invoked once all FASTA files have been read; the default implementation does nothing.
     *
     * @param contigsPerTaxid the trie tracking how many contigs were seen per taxid
     */
    protected void afterReadFastas(AbstractRefSeqFastaReader.StringLong2DigitTrie contigsPerTaxid) {
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
            // Copied before iterating: a consumer still dying adds to the list while we walk it.
            List<Throwable> throwables = new ArrayList<>(bundle.getThrowableList());
            for (Throwable t : throwables) {
                if (getLogger().isErrorEnabled()) {
                    getLogger().error("Error in consumer thread: ", t);
                }
            }
            bundle.clearThrowableList();
            // The first one as the cause, as AbstractFastqReader does it: a caller that catches this
            // and prints only the message would otherwise be left with nothing to act on.
            throw throwables.isEmpty() ? new RuntimeException("Error(s) in consumer thread(s).")
                    : new RuntimeException("Error(s) in consumer thread(s).", throwables.get(0));
        }
    }

    /**
     * Creates a consumer {@link Runnable} that takes files from the queue and reads them with its
     * own FASTA reader.
     *
     * @param i               the index of the consumer thread
     * @param blockingQueue   the queue supplying files to read
     * @param contigsPerTaxid the shared trie tracking how many contigs were seen per taxid
     * @return the consumer runnable
     */
    protected Runnable createFastaReaderRunnable(int i,
                                                 BlockingQueue<FileAndNode> blockingQueue,
                                                 AbstractRefSeqFastaReader.StringLong2DigitTrie contigsPerTaxid) {
        AbstractRefSeqFastaReader fastaReader = createFastaReader(contigsPerTaxid);
        return new Runnable() {
            @Override
            public void run() {
                // A pass is ended by interrupting the consumers, and a thread that was not waiting at
                // that moment carries the interrupt back into the pool with it. Clearing it here keeps
                // it from striking the next pass, which would take it for a failure of its own.
                Thread.interrupted();
                while (!dump) {
                    try {
                        try {
                            FileAndNode fileAndNode = blockingQueue.take();
                            if (fileAndNode == END_OF_PASS) {
                                // Returned, not broken out of: the enclosing finally would otherwise
                                // count down for a file that was never counted up.
                                return;
                            }
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
     * Creates the FASTA reader that processes each contig; called once per reader thread. The shared
     * {@code contigsPerTaxid} trie tracks how many contigs have been seen per taxid.
     *
     * @param contigsPerTaxid the shared trie tracking how many contigs were seen per taxid
     * @return the FASTA reader that processes each contig
     */
    protected abstract AbstractRefSeqFastaReader createFastaReader(AbstractRefSeqFastaReader.StringLong2DigitTrie contigsPerTaxid);

    /**
     * In addition to discarding the result, aborts a pass that may still be reading.
     */
    public void dump() {
        super.dump();
        cleanUpThreads();
    }

    /**
     * Aborts the pass that is reading, if one is: signals its consumer threads to stop and interrupts
     * any that are blocked.
     * <p>
     * For aborting only. A pass that has read everything ends its own consumers - see
     * {@link #endOfPass(BlockingQueue)} - and calling this at the end of one instead, as every
     * subclass used to do from its {@code finally}, leaves the flag set behind: nothing clears it on
     * the way out, so the goal reads nothing at all the next time it is made, and reports success for
     * it. Note also that {@link ExecutionContext#interruptAll()} interrupts the calling thread along
     * with the consumers, which a caller that goes on to read again has to survive.
     */
    protected void cleanUpThreads() {
        dump = true;
        bundle.interruptAll();
    }

    /**
     * Clears the abort flag so that this goal can read again.
     * <p>
     * {@link #readFastas()} calls this itself, at the start of every pass, so that a goal aborted
     * through {@link #dump()} is not thereby finished for good. Nothing else has to call it: a pass
     * that ends normally never sets the flag in the first place.
     */
    protected void readyForAnotherPass() {
        dump = false;
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
