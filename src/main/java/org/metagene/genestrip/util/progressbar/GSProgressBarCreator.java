package org.metagene.genestrip.util.progressbar;

import me.tongfei.progressbar.*;
import org.apache.commons.logging.Log;

import java.time.Duration;
import java.util.Optional;

public class GSProgressBarCreator {
    private GSProgressBarCreator() {}

    public static ProgressBar newGSProgressBar(String task, GSProgressUpdate progressUpdate, Log log) {
        return newGSProgressBar(task, 0,1000, " bytes", progressUpdate, log);
    }

    public static ProgressBar newGSProgressBar(String task, int max, String unitName, Log log) {
        return newGSProgressBar(task, max, 1000, unitName, null, log);
    }

    public static ProgressBar newGSProgressBar(String task, String unitName, GSProgressUpdate progressUpdate, Log log) {
        return newGSProgressBar(task, 0, 1000, unitName, progressUpdate, log);
    }

    public static ProgressBar newGSProgressBar(String task, int max, int updateIntervalMillis, String unitName, GSProgressUpdate progressUpdate, Log log) {
        GSProgressBarRenderer renderer = new GSProgressBarRenderer(unitName, progressUpdate);
        ProgressBarBuilder progressBarBuilder = new ProgressBarBuilder().setInitialMax(max).
                setTaskName(task).setUpdateIntervalMillis(updateIntervalMillis).
                setUnit(unitName, 1).setRenderer(renderer).setMaxRenderedLength(100).continuousUpdate();
        if (log != null && log.isInfoEnabled()) {
            progressBarBuilder.setConsumer(new DelegatingProgressBarConsumer(log::info));
        }
        ProgressBar progressBar = progressBarBuilder.build();
        renderer.setProgressBar(progressBar);
        return progressBar;
    }

    private static class GSProgressBarRenderer extends DefaultProgressBarRenderer {
        private final GSProgressUpdate progressUpdate;
        private ProgressBar progressBar;

        protected GSProgressBarRenderer(String unitName, GSProgressUpdate progressUpdateProvider) {
            super(ProgressBarStyle.ASCII, unitName, 1, false, null, null, true, GSProgressBarRenderer::linearEta);
            this.progressUpdate = progressUpdateProvider;
        }

        public void setProgressBar(ProgressBar progressBar) {
            this.progressBar = progressBar;
        }

        @Override
        public String render(ProgressState progress, int maxLength) {
            if (progressUpdate != null) {
                progressBar.maxHint(progressUpdate.max());
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
