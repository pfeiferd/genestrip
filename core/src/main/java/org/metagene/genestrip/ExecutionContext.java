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
package org.metagene.genestrip;

import java.util.List;

import org.metagene.genestrip.io.StreamingResource;

/**
 * Manages the worker threads used for parallel data processing and collects any throwables raised
 * within them. Implementations also receive progress callbacks and control interruption of running
 * work.
 */
public interface ExecutionContext {
	/**
	 * Returns the number of worker threads configured for parallel processing.
	 *
	 * @return the configured number of worker threads
	 */
	public int getThreads();

	/**
	 * Executes the given task, either on a worker thread or inline in the calling thread depending
	 * on the configured number of threads.
	 *
	 * @param runnable the task to execute
	 */
	public void execute(Runnable runnable);

	/**
	 * Returns the list of throwables raised by worker threads so far.
	 *
	 * @return the list of collected throwables
	 */
	public List<Throwable> getThrowableList();

	/**
	 * Returns whether any throwables have been raised by worker threads.
	 *
	 * @return {@code true} if at least one throwable has been collected
	 */
	public boolean hasThrowables();

	/**
	 * Clears the list of collected throwables.
	 */
	public void clearThrowableList();

	/**
	 * Interrupts all worker threads (and the main thread, if any) to abort ongoing processing.
	 */
	public void interruptAll();

	/**
	 * Reports overall progress across all inputs currently being processed.
	 *
	 * @param coveredBytes the number of bytes processed so far
	 * @param estTotalBytes the estimated total number of bytes to process
	 * @param elapsedTimeMS the elapsed processing time in milliseconds
	 * @param estTotalTimeMS the estimated total processing time in milliseconds
	 * @param ratio the fraction of work completed
	 * @param totalProcessedReads the total number of reads processed so far
	 * @param index the index of the input currently being processed
	 * @param totalCount the total number of inputs to process
	 */
	public void setTotalProgress(long coveredBytes, long estTotalBytes, long elapsedTimeMS, long estTotalTimeMS,
			double ratio, long totalProcessedReads, int index, int totalCount);

	/**
	 * Reports progress for a single fastq/fasta resource being processed.
	 *
	 * @param fastq the resource currently being processed
	 * @param coveredBytes the number of bytes of the resource processed so far
	 * @param estTotalBytes the estimated total number of bytes in the resource
	 * @param elapsedTimeMS the elapsed processing time in milliseconds
	 * @param estTotalTimeMS the estimated total processing time in milliseconds
	 * @param ratio the fraction of the resource completed
	 * @param processedReads the number of reads processed from the resource so far
	 */
	public void setFastqProgress(StreamingResource fastq, long coveredBytes, long estTotalBytes,
			long elapsedTimeMS, long estTotalTimeMS, double ratio, long processedReads);

	/**
	 * Returns whether this context wants progress callbacks.
	 *
	 * @return whether this context wants progress callbacks (e.g. to drive a progress bar)
	 */
	public boolean isRequiresProgress();

	/**
	 * Returns the interval, in processed units, between progress log updates.
	 *
	 * @return the log update cycle length
	 */
	public long getLogUpdateCycle();
	
	/**
	 * Releases resources held by this context, in particular shutting down the thread pool.
	 */
	public void dump();
}