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
package org.metagene.genestrip.fastq;

import java.io.IOException;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.Locale;

import org.metagene.genestrip.ExecutionContext;
import org.metagene.genestrip.io.StreamingResource;
import org.metagene.genestrip.io.StreamingResourceStream;
import me.tongfei.progressbar.*;
import org.metagene.genestrip.util.GSLogFactory;
import org.metagene.genestrip.util.progressbar.GSProgressBarCreator;


/**
 * An {@link AbstractFastqReader} that reads through a stream of FASTQ/FASTA resources one after
 * another while logging and reporting per-file and overall progress.
 */
public abstract class AbstractLoggingFastqStreamer extends AbstractFastqReader {
    /** Type hint identifying a resource that should be parsed as FASTA. */
    public static final String FASTA_TYPE_HINT = "fasta";
    /** Type hint identifying a resource that should be parsed as FASTQ. */
    public static final String FASTQ_TYPE_HINT = "fastq";

    private static final DecimalFormat DF = new DecimalFormat("0.0", new DecimalFormatSymbols(Locale.US));

    private StreamingResource.StreamAccess byteCountAccess;
    private int coveredCounter;
    private int totalCount;
    private StreamingResource currentFastq;
    private long fastqsFileSize;
    private long fastqStartTime;
    private long fastqFileSize;
    private long coveredFilesSize;
    private long startTime;
    /** The number of reads indexed (written) so far. */
    protected long indexedC;
    /** The total number of reads processed across all files. */
    protected long totalReads;
    /** The total number of k-mers processed across all files. */
    protected long totalKMers;
    /** The total number of base pairs processed across all files. */
    protected long totalBPs;
    private final long logUpdateCycle;

    /**
     * Creates the streamer with the given matching parameters and execution context.
     *
     * @param k the k-mer length
     * @param initialReadSize the initial read buffer size in bytes
     * @param maxQueueSize the maximum size of the consumer thread queue
     * @param bundle the execution context providing threads and progress reporting
     * @param withProbs whether match probabilities are computed
     * @param config additional configuration passed to the underlying reader
     */
    public AbstractLoggingFastqStreamer(int k, int initialReadSize, int maxQueueSize, ExecutionContext bundle, boolean withProbs,
                                        Object... config) {
        super(k, initialReadSize, maxQueueSize, bundle, withProbs, config);

        this.logUpdateCycle = bundle.getLogUpdateCycle();
    }

    /**
     * Reads and processes every FASTQ/FASTA resource in the given stream, logging and reporting
     * progress and accumulating total read, k-mer and base-pair counts.
     *
     * @param fastqs the stream of FASTQ/FASTA resources to process
     * @throws IOException if a resource cannot be opened or read
     */
    public void processFastqStreams(StreamingResourceStream fastqs) throws IOException {
        if (logger.isInfoEnabled()) {
            logger.info("Number of consumer threads: " + bundle.getThreads());
        }

        startTime = System.currentTimeMillis();
        totalCount = fastqs.size();
        fastqsFileSize = fastqs.getTotalByteSize();
        coveredFilesSize = 0;
        coveredCounter = 0;
        totalReads = 0;

        for (StreamingResource fastq : fastqs) {
            currentFastq = fastq;
            try (StreamingResource.StreamAccess lbyteCountAccess = fastq.openStream()) {
                byteCountAccess = lbyteCountAccess;
                fastqFileSize = byteCountAccess.getSize();
                if (fastqFileSize != -1 && fastqsFileSize == -1) {
                    if (coveredCounter == totalCount - 1) {
                        fastqsFileSize = coveredFilesSize + fastqFileSize;
                    }
                }
                fastqStartTime = System.currentTimeMillis();
                logFastqStart();
                try (ProgressBar pb = isProgressBar() ? GSProgressBarCreator.newGSProgressBar(getProgressBarTaskName(), byteCountAccess, null) : null) {
                    readFastq(byteCountAccess.getInputStream(), isFastaStream(fastq));
                }
                logFastqDone();
                totalReads += reads;
                totalKMers += kMers;
                totalBPs += readBPs;
                coveredFilesSize += byteCountAccess.getBytesRead();
            }
            coveredCounter++;
        }
        allDone();
    }

    /**
     * Returns whether the given resource should be parsed as FASTA (based on its type hint).
     *
     * @param fastq the resource to inspect
     * @return whether the given resource should be parsed as FASTA
     */
    protected boolean isFastaStream(StreamingResource fastq) {
        return FASTA_TYPE_HINT.equals(fastq.getTypeHint());
    }

    /**
     * Returns whether a progress bar should be shown while processing.
     *
     * @return whether a progress bar should be shown
     */
    protected boolean isProgressBar() {
        return true;
    }

    /**
     * Returns the task name displayed on the progress bar.
     *
     * @return the progress bar task name
     */
    protected String getProgressBarTaskName() {
        return ((GSLogFactory.GSLog) logger).getName();
    }

    @Override
    protected void start() throws IOException {
        indexedC = 0;
        doUpdateProgress();
    }

    @Override
    protected void done() throws IOException {
        if (bundle.isRequiresProgress()) {
            long bytesCovered = byteCountAccess.getBytesRead();
            long totalTime = System.currentTimeMillis() - fastqStartTime;
            bundle.setFastqProgress(currentFastq, bytesCovered, fastqFileSize, totalTime, totalTime, 1, reads);
        }
    }

    @Override
    protected void updateWriteStats() {
        indexedC++;
    }

    @Override
    protected void updateProgress() {
        if (logUpdateCycle > 0 && reads % logUpdateCycle == 0) {
            doUpdateProgress();
        }
    }

    /**
     * Computes and reports progress for the current file and for the overall run to the execution
     * context.
     */
    protected void doUpdateProgress() {
        if (bundle.isRequiresProgress()) {
            long bytesCovered = byteCountAccess.getBytesRead();
            double ratio = fastqFileSize == -1 ? -1 : bytesCovered / (double) fastqFileSize;
            long stopTime = System.currentTimeMillis();
            long diff = (stopTime - fastqStartTime);
            double totalTime = ratio == -1 ? -1 : diff / ratio;
            bundle.setFastqProgress(currentFastq, bytesCovered, fastqFileSize, diff, (long) totalTime, ratio,
                    reads);

            bytesCovered = (coveredFilesSize + byteCountAccess.getBytesRead());
            ratio = fastqsFileSize == -1 ? -1 : bytesCovered / (double) fastqsFileSize;
            diff = (stopTime - startTime);
            totalTime = ratio == -1 ? -1 : diff / ratio;
            bundle.setTotalProgress(bytesCovered, fastqsFileSize, diff, (long) totalTime, ratio, totalReads + reads,
                    coveredCounter, totalCount);
        }
    }

    /**
     * Returns the number of reads between progress log updates.
     *
     * @return the log update cycle in reads
     */
    public long getLogUpdateCycle() {
        return logUpdateCycle;
    }

    /**
     * Logs the start of processing for the current file.
     */
    protected void logFastqStart() {
        if (logger.isInfoEnabled()) {
            logger.info("Processing fastq file (" + (coveredCounter + 1) + "/" + totalCount + "): " + currentFastq);
        }
    }

    /**
     * Logs completion statistics for the current file.
     */
    protected void logFastqDone() {
        if (logger.isInfoEnabled()) {
            long bytesCovered = byteCountAccess.getBytesRead();
            long totalTime = System.currentTimeMillis() - fastqStartTime;
            double totalHours = ((double) totalTime) / 1000 / 60 / 60;
            logger.info("Done with fastq: " + currentFastq);
            logger.info("Hours: " + DF.format(totalHours));
            logger.info("Bytes: " + bytesCovered);
            logger.info("Reads: " + reads);
        }
    }

    /**
     * Logs and reports the final totals after all files have been processed.
     */
    protected void allDone() {
        long totalTime = (System.currentTimeMillis() - startTime);
        if (bundle.isRequiresProgress()) {
            bundle.setTotalProgress(coveredFilesSize, coveredFilesSize, totalTime, totalTime, 1, totalReads,
                    coveredCounter, totalCount);
        }
        if (logger.isInfoEnabled()) {
            double totalHours = ((double) totalTime) / 1000 / 60 / 60;
            logger.info("All done with fastqs.");
            logger.info("Total hours: " + DF.format(totalHours));
            logger.info("Total bytes: " + coveredFilesSize);
            logger.info("Total reads: " + totalReads);
        }
    }
}
