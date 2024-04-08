package org.metagene.genestrip.util;

import java.lang.Thread.UncaughtExceptionHandler;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.ThreadFactory;

public class DefaultExecutorServiceBundle implements ExecutorServiceBundle {
	private final int consumers;
	private final List<Throwable> throwablesInThreads;
	private final List<Throwable> immutableThrowablesInThreads;
	private final ExecutorService executorService;
	private final List<Thread> executorThreads;
	private final long logUpdateCycle;
	
	private long elapsedTimeMS;
	private long estTotalTimeMS;
	private long coveredBytes;
	private long estTotalBytes;
	private double progressRatio;

	public DefaultExecutorServiceBundle(int consumers, long logUpdateCycle) {
		this.consumers = consumers;
		this.throwablesInThreads = Collections.synchronizedList(new ArrayList<Throwable>());
		this.logUpdateCycle = logUpdateCycle;
		immutableThrowablesInThreads = Collections.unmodifiableList(throwablesInThreads);
		executorThreads = new ArrayList<Thread>();
		this.executorService = consumers > 0 ? Executors.newFixedThreadPool(consumers, createThreadFactory()) : null;
	}
	
	@Override
	public long getLogUpdateCycle() {
		return logUpdateCycle;
	}
	
	@Override
	public void setProgress(long coveredBytes, long estTotalBytes, long elapsedTimeMS, long estTotalTimeMS, 
			double ratio) {
		this.coveredBytes = coveredBytes;
		this.estTotalBytes = estTotalBytes;
		this.elapsedTimeMS = elapsedTimeMS;
		this.estTotalTimeMS = estTotalTimeMS;
		this.progressRatio = ratio;
	}

	
	@Override
	public long getElapsedTimeMS() {
		return elapsedTimeMS;
	}
	
	@Override
	public long getEstTotalTimeMS() {
		return estTotalTimeMS;
	}
	
	@Override
	public double getProgressRatio() {
		return progressRatio;
	}
	
	@Override
	public long getCoveredBytes() {
		return coveredBytes;
	}
	
	@Override
	public long getEstTotalBytes() {
		return estTotalBytes;
	}
	
	@Override
	public boolean isRequiresProgress() {
		return false;
	}
	
	@Override
	public void interruptAll() {
		for (Thread t : executorThreads) {
			t.interrupt();
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
	public void clearThrowableList() {
		throwablesInThreads.clear();
	}

	protected ThreadFactory createThreadFactory() {
		return new ThreadFactory() {
			private int newCounter = 0;

			@Override
			public Thread newThread(Runnable r) {
				Thread t = new Thread(r);
				t.setName(getThreadBaseName() + " #" + newCounter++);
				t.setUncaughtExceptionHandler(new UncaughtExceptionHandler() {
					@Override
					public void uncaughtException(Thread t, Throwable e) {
						throwablesInThreads.add(e);
					}
				});
				executorThreads.add(t);
				return t;
			}
		};
	}

	protected String getThreadBaseName() {
		return "Executor Thread";
	}
	
	@Override
	public void dump() {
		executorService.shutdown();
	}
}