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


public abstract class AbstractLoggingFastqStreamer extends AbstractFastqReader {
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
    protected long indexedC;
    protected long totalReads;
    protected long totalKMers;
    protected long totalBPs;
    private final long logUpdateCycle;

    public AbstractLoggingFastqStreamer(int k, int initialReadSize, int maxQueueSize, ExecutionContext bundle, boolean withProbs,
                                        Object... config) {
        super(k, initialReadSize, maxQueueSize, bundle, withProbs, config);

        this.logUpdateCycle = bundle.getLogUpdateCycle();
    }

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
                    readFastq(byteCountAccess.getInputStream());
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

    protected boolean isProgressBar() {
        return true;
    }

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

    public long getLogUpdateCycle() {
        return logUpdateCycle;
    }

    protected void logFastqStart() {
        if (logger.isInfoEnabled()) {
            logger.info("Processing fastq file (" + (coveredCounter + 1) + "/" + totalCount + "): " + currentFastq);
        }
    }

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
