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
package org.metagene.genestrip.goals.genbank;

import java.io.File;
import java.util.Arrays;
import java.util.List;

import org.metagene.genestrip.GSGoalKey;
import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.genbank.AssemblySummaryReader;
import org.metagene.genestrip.goals.GSFileDownloadGoal;
import org.metagene.genestrip.make.Goal;

public class AssemblyFileDownloadGoal<P extends GSProject> extends GSFileDownloadGoal<P> {
	public static final String FTP_DIR_GENBANK = "/genomes/genbank";
	public static final String FTP_DIR_REFSEQ = "/genomes/refseq";

	private final List<File> files;

	@SafeVarargs
	public AssemblyFileDownloadGoal(P project, Goal<P>... deps) {
		super(project, GSGoalKey.ASSEMBLYDOWNLOAD, deps);
		files = Arrays.asList(
				new File(project.getCommon().getGenbankDir(), AssemblySummaryReader.ASSEMLY_SUM_GENBANK),
				new File(project.getCommon().getRefSeqDir(), AssemblySummaryReader.ASSEMLY_SUM_REFSEQ));
	}

	@Override
	public List<File> getFiles() {
		return files;
	}

	@Override
	protected String getFTPDir(File file) {
		return file.getName().endsWith(AssemblySummaryReader.ASSEMLY_SUM_GENBANK) ? FTP_DIR_GENBANK : FTP_DIR_REFSEQ;
	}
}
