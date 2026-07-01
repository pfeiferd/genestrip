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

import java.io.File;
import java.util.Collections;
import java.util.List;

import org.metagene.genestrip.GSGoalKey;
import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.make.Goal;

/**
 * Download goal for the RefSeq {@code RELEASE_NUMBER} file, which identifies the current RefSeq
 * release.
 *
 * @param <P> the project type
 */
public class RefSeqRNumDownloadGoal<P extends GSProject> extends RefSeqDownloadGoal<P> {
	/** The name of the RefSeq release-number file. */
	public static final String RELEASE_NUMBER_FILE_NAME = "RELEASE_NUMBER";

	private final List<File> files;

	/**
	 * Creates the download goal for the RefSeq release-number file.
	 *
	 * @param project the project this goal belongs to.
	 * @param deps the goals this goal depends on.
	 */
	@SafeVarargs
	public RefSeqRNumDownloadGoal(P project, Goal<P>... deps) {
		super(project, GSGoalKey.REFSEQRELEASE, deps);
		
		files = Collections.singletonList(new File(project.getCommon().getRefSeqDir(), RELEASE_NUMBER_FILE_NAME));
	}

	@Override
	protected String getFTPDir(File file) {
		return RELEASE_FOLDER;
	}
	
	@Override
	public List<File> getFiles() {
		return files;
	}
}
