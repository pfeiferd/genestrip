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
package org.metagene.genestrip.goals;

import me.tongfei.progressbar.DelegatingProgressBarConsumer;
import me.tongfei.progressbar.ProgressBar;
import me.tongfei.progressbar.ProgressBarBuilder;
import me.tongfei.progressbar.ProgressBarStyle;
import org.metagene.genestrip.GSConfigKey;
import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.make.FileDownloadGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.GoalKey;
import org.metagene.genestrip.util.progressbar.GSProgressBarCreator;

public abstract class GSFileDownloadGoal extends FileDownloadGoal<GSProject> {
	@SafeVarargs
	public GSFileDownloadGoal(GSProject project, GoalKey key, Goal<GSProject>... deps) {
		super(project, key, deps);
	}

	@Override
	protected String getHttpBaseURL() {
		return stringConfigValue(GSConfigKey.HTTP_BASE_URL);
	}

	@Override
	protected String getFTPBaseURL() {
		return stringConfigValue(GSConfigKey.FTP_BASE_URL);
	}

	@Override
	protected boolean isUseHttp() {
		return booleanConfigValue(GSConfigKey.USE_HTTP);
	}

	@Override
	protected boolean isIgnoreMissingFiles() {
		return booleanConfigValue((GSConfigKey.IGNORE_MISSING_FASTAS));
	}

	@Override
	protected int getMaxDownloadTries() {
		return intConfigValue(GSConfigKey.MAX_DOWNLOAD_TRIES);
	}

	@Override
	protected ProgressBar createProgressBar(int max) {
		return booleanConfigValue(GSConfigKey.PROGRESS_BAR) ?
				GSProgressBarCreator.newGSProgressBar(getKey().getName(), max, 60000, " files", null, getLogger()) :
				null;
	}
}
