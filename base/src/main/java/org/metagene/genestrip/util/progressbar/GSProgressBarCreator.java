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
package org.metagene.genestrip.util.progressbar;

import me.tongfei.progressbar.*;
import org.apache.commons.logging.Log;

import java.time.Duration;
import java.util.Optional;

/**
 * Factory for tongfei {@link me.tongfei.progressbar.ProgressBar} instances configured for Genestrip:
 * ASCII style, a custom renderer that pulls live current/max values from a {@link GSProgressUpdate},
 * a configurable update interval and unit name, and optional rendering of the bar to a {@link Log}
 * at debug level.
 */
public class GSProgressBarCreator {
    private static int globalUpdateInterval = 1000;

    /**
     * Sets the default refresh interval (in milliseconds) used for newly created progress bars.
     *
     * @param globalUpdateInterval the refresh interval in milliseconds
     */
    public static void setGlobalUpdateInterval(int globalUpdateInterval) {
        GSProgressBarCreator.globalUpdateInterval = globalUpdateInterval;
    }

    /**
     * Returns the default refresh interval (in milliseconds) used for newly created progress bars.
     *
     * @return the refresh interval in milliseconds
     */
    public static int getGlobalUpdateInterval() {
        return globalUpdateInterval;
    }

    private GSProgressBarCreator() {}

    /**
     * Creates a byte-counting progress bar for the given task that reads its current and maximum
     * values from the given progress update.
     *
     * @param task           the task name shown on the bar
     * @param progressUpdate a live source of current/max values
     * @param log            a log to render the bar to at debug level, or {@code null} for the console
     * @return the created progress bar
     */
    public static ProgressBar newGSProgressBar(String task, GSProgressUpdate progressUpdate, Log log) {
        return newGSProgressBar(task, 0,globalUpdateInterval, " bytes", progressUpdate, log, false);
    }

    /**
     * Creates a byte-counting progress bar as {@link #newGSProgressBar(String, GSProgressUpdate, Log)}.
     * If {@code doneWorkaround} is {@code true}, the bar jumps to its maximum once the task is no
     * longer alive.
     *
     * @param task           the task name shown on the bar
     * @param progressUpdate a live source of current/max values
     * @param log            a log to render the bar to at debug level, or {@code null} for the console
     * @param doneWorkaround whether the bar jumps to its maximum once the task is no longer alive
     * @return the created progress bar
     */
    public static ProgressBar newGSProgressBar(String task, GSProgressUpdate progressUpdate, Log log, boolean doneWorkaround) {
        return newGSProgressBar(task, 0,globalUpdateInterval, " bytes", progressUpdate, log, doneWorkaround);
    }

    /**
     * Creates a progress bar for the given task with a fixed maximum and the given unit name, stepped
     * by the caller.
     *
     * @param task     the task name shown on the bar
     * @param max      the fixed maximum value
     * @param unitName the unit label for counted amounts
     * @param log      a log to render the bar to at debug level, or {@code null} for the console
     * @return the created progress bar
     */
    public static ProgressBar newGSProgressBar(String task, long max, String unitName, Log log) {
        return newGSProgressBar(task, max, globalUpdateInterval, unitName, null, log, false);
    }

    /**
     * Creates a progress bar for the given task with the given unit name that reads its progress from
     * the given progress update.
     *
     * @param task           the task name shown on the bar
     * @param unitName       the unit label for counted amounts
     * @param progressUpdate a live source of current/max values
     * @param log            a log to render the bar to at debug level, or {@code null} for the console
     * @return the created progress bar
     */
    public static ProgressBar newGSProgressBar(String task, String unitName, GSProgressUpdate progressUpdate, Log log) {
        return newGSProgressBar(task, 0, globalUpdateInterval, unitName, progressUpdate, log, false);
    }

    /**
     * Creates a progress bar with full control over its parameters.
     *
     * @param task                 the task name shown on the bar
     * @param max                  the initial maximum value
     * @param updateIntervalMillis how often the bar is refreshed, in milliseconds
     * @param unitName             the unit label for counted amounts
     * @param progressUpdate       a live source of current/max values, or {@code null} to step the bar
     *                             manually
     * @param log                  a log to render the bar to at debug level, or {@code null} for the
     *                             console
     * @param doneWorkaround       whether the bar jumps to its maximum once the task is no longer alive
     * @return the created progress bar
     */
    public static ProgressBar newGSProgressBar(String task, long max, int updateIntervalMillis, String unitName, GSProgressUpdate progressUpdate, Log log, boolean doneWorkaround) {
        GSProgressBarRenderer renderer = new GSProgressBarRenderer(unitName, progressUpdate, doneWorkaround);
        ProgressBarBuilder progressBarBuilder = new ProgressBarBuilder().setInitialMax(max).
                setTaskName(task).setUpdateIntervalMillis(updateIntervalMillis).
                setUnit(unitName, 1).setRenderer(renderer).setMaxRenderedLength(100).continuousUpdate();
        if (log != null && log.isDebugEnabled()) {
            progressBarBuilder.setConsumer(new DelegatingProgressBarConsumer(log::debug));
        }
        ProgressBar progressBar = progressBarBuilder.build();
        renderer.setProgressBar(progressBar);
        return progressBar;
    }

    private static class GSProgressBarRenderer extends DefaultProgressBarRenderer {
        private final GSProgressUpdate progressUpdate;
        private final boolean doneWorkaround;
        private ProgressBar progressBar;

        protected GSProgressBarRenderer(String unitName, GSProgressUpdate progressUpdateProvider, boolean doneWorkaround) {
            super(ProgressBarStyle.ASCII, unitName, 1, false, null, null, true, GSProgressBarRenderer::linearEta);
            this.progressUpdate = progressUpdateProvider;
            this.doneWorkaround = doneWorkaround;
        }

        public void setProgressBar(ProgressBar progressBar) {
            this.progressBar = progressBar;
        }

        @Override
        public String render(ProgressState progress, int maxLength) {
            if (progressUpdate != null) {
                progressBar.maxHint(progressUpdate.max());
                progressBar.stepTo(!progress.isAlive() && doneWorkaround ? progressUpdate.max() : progressUpdate.current());
            }
            return super.render(progress, maxLength);
        }

        static Optional<Duration> linearEta(ProgressState progress) {
            if (progress.getMax() <= 0 || progress.isIndefinite()) return Optional.empty();
            else if (progress.getCurrent() - progress.getStart() == 0) return Optional.empty();
            else return Optional.of(
                        progress.getElapsedAfterStart()
                                .dividedBy(progress.getCurrent() - progress.getStart())
                                .multipliedBy(progress.getMax() - progress.getCurrent())
                );
        }
    }
}
