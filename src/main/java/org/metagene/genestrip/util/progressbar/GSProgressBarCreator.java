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

public class GSProgressBarCreator {
    private static int globalUpdateInterval = 1000;

    public static void setGlobalUpdateInterval(int globalUpdateInterval) {
        GSProgressBarCreator.globalUpdateInterval = globalUpdateInterval;
    }

    public static int getGlobalUpdateInterval() {
        return globalUpdateInterval;
    }

    private GSProgressBarCreator() {}

    public static ProgressBar newGSProgressBar(String task, GSProgressUpdate progressUpdate, Log log) {
        return newGSProgressBar(task, 0,globalUpdateInterval, " bytes", progressUpdate, log, false);
    }

    public static ProgressBar newGSProgressBar(String task, GSProgressUpdate progressUpdate, Log log, boolean doneWorkaround) {
        return newGSProgressBar(task, 0,globalUpdateInterval, " bytes", progressUpdate, log, doneWorkaround);
    }

    public static ProgressBar newGSProgressBar(String task, int max, String unitName, Log log) {
        return newGSProgressBar(task, max, globalUpdateInterval, unitName, null, log, false);
    }

    public static ProgressBar newGSProgressBar(String task, String unitName, GSProgressUpdate progressUpdate, Log log) {
        return newGSProgressBar(task, 0, globalUpdateInterval, unitName, progressUpdate, log, false);
    }

    public static ProgressBar newGSProgressBar(String task, int max, int updateIntervalMillis, String unitName, GSProgressUpdate progressUpdate, Log log, boolean doneWorkaround) {
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
                progressBar.maxHint(!progress.isAlive() && doneWorkaround ? progressUpdate.current() : progressUpdate.max());
                progressBar.stepTo(progressUpdate.current());
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
