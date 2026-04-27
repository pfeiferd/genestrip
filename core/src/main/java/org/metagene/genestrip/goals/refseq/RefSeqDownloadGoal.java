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
package org.metagene.genestrip.goals.refseq;

import org.metagene.genestrip.GSConfigKey;
import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.goals.GSFileDownloadGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.GoalKey;

public abstract class RefSeqDownloadGoal<P extends GSProject> extends GSFileDownloadGoal<P> {
	public static final String RELEASE_FOLDER = "/release";

	@SafeVarargs
	public RefSeqDownloadGoal(P project, GoalKey key, Goal<P>... deps) {
		super(project, key, deps);
	}
	
	@Override
	protected String getHttpBaseURL() {
		return stringConfigValue(GSConfigKey.REF_SEQ_HTTP_BASE_URL);
	}
	
	@Override
	protected String getFTPBaseURL() {
		return stringConfigValue(GSConfigKey.REF_SEQ_FTP_BASE_URL);
	}
}
