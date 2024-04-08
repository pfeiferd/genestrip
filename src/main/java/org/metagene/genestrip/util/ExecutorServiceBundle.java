package org.metagene.genestrip.util;

import java.util.List;

public interface ExecutorServiceBundle {
	public int getThreads();

	public void execute(Runnable runnable);

	public List<Throwable> getThrowableList();

	public void clearThrowableList();

	public void interruptAll();

	public void setProgress(long coveredBytes, long estTotalBytes, long elapsedTimeMS, long estTotalTimeMS,
			double ratio);

	public boolean isRequiresProgress();

	public long getElapsedTimeMS();

	public long getEstTotalTimeMS();

	public double getProgressRatio();

	public long getCoveredBytes();

	public long getEstTotalBytes();

	public long getLogUpdateCycle();
	
	public void dump();
}