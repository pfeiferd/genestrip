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

import java.lang.Thread.UncaughtExceptionHandler;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.ThreadFactory;

import org.metagene.genestrip.io.StreamingResource;

public class DefaultExecutionContext implements ExecutionContext {
	private final Thread mainThread;
	private final String threadBaseName;
	private final int consumers;
	private boolean hasThrowables;
	private final List<Throwable> throwablesInThreads;
	private final List<Throwable> immutableThrowablesInThreads;
	private final ExecutorService executorService;
	private final List<Thread> executorThreads;
	private final long logUpdateCycle;

	public DefaultExecutionContext(Thread mainTread, int consumers, long logUpdateCycle) {
		this(mainTread, consumers, logUpdateCycle, "Executor Thread");
	}
	
	public DefaultExecutionContext(Thread mainTread, int consumers, long logUpdateCycle, String threadBaseName) {
		this.mainThread = mainTread;
		this.threadBaseName = threadBaseName;
		this.consumers = consumers < 0 ? (Runtime.getRuntime().availableProcessors() - 1) : consumers;
		this.throwablesInThreads = Collections.synchronizedList(new ArrayList<Throwable>());
		this.logUpdateCycle = logUpdateCycle;
		immutableThrowablesInThreads = Collections.unmodifiableList(throwablesInThreads);
		executorThreads = new ArrayList<Thread>();
		this.executorService = this.consumers > 0 ? Executors.newFixedThreadPool(this.consumers, createThreadFactory()) : null;
	}
	
	@Override
	public long getLogUpdateCycle() {
		return logUpdateCycle;
	}
	
	@Override
	public void setTotalProgress(long coveredBytes, long estTotalBytes, long elapsedTimeMS, long estTotalTimeMS,
			double ratio, long totalProcessedReads, int index, int totalCount) {
	}

	@Override
	public void setFastqProgress(StreamingResource fastq, long coveredBytes, long estTotalBytes,
			long elapsedTimeMS, long estTotalTimeMS, double ratio, long processedReads) {		
	}
	
	@Override
	public boolean isRequiresProgress() {
		return false;
	}
		
	@Override
	public void interruptAll() {
		// Create a copy to avoid comodifaction problems in case of error situations:
		List<Thread> executorThreads = new ArrayList<>(this.executorThreads);
		for (Thread t : executorThreads) {
			t.interrupt();
		}
		if (mainThread != null) {
			mainThread.interrupt();
		}
	}

	@Override
	public int getThreads() {
		return consumers;
	}

	@Override
	public void execute(Runnable runnable) {
		if (executorService != null) {
			executorService.execute(runnable);
		}
		else {
			runnable.run();
		}
	}

	@Override
	public List<Throwable> getThrowableList() {
		return immutableThrowablesInThreads;
	}

	@Override
	public final boolean hasThrowables() {
		return hasThrowables;
	}

	@Override
	public void clearThrowableList() {
		throwablesInThreads.clear();
		hasThrowables = false;
	}

	protected ThreadFactory createThreadFactory() {
		return new ThreadFactory() {
			private int newCounter = 0;

			@Override
			public Thread newThread(Runnable r) {
				Thread t = new Thread(r);
				t.setName(threadBaseName + " #" + newCounter++);
				t.setUncaughtExceptionHandler(new UncaughtExceptionHandler() {
					@Override
					public void uncaughtException(Thread t, Throwable e) {
						throwablesInThreads.add(e);
						hasThrowables = true;
						interruptAll();
					}
				});
				executorThreads.add(t);
				return t;
			}
		};
	}
	
	@Override
	public void dump() {
		if (executorService != null) {
			executorService.shutdown();
		}
	}
}