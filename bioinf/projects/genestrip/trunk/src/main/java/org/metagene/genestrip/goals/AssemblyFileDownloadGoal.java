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

import java.io.File;
import java.util.Collections;
import java.util.List;

import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.make.FileDownloadGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.tax.AssemblySummaryReader;

public class AssemblyFileDownloadGoal extends FileDownloadGoal<GSProject> {
	public static final String FTP_DIR_GENBANK = "/genomes/genbank";
	public static final String FTP_DIR_REFSEQ = "/genomes/refseq";

	private final List<File> files;

	@SafeVarargs
	public AssemblyFileDownloadGoal(GSProject project, String name, Goal<GSProject>... deps) {
		super(project, name, deps);
		files = Collections.singletonList(new File(project.getConfig().getCommonDir(),
				project.getConfig().isUseGenBank() ? AssemblySummaryReader.ASSEMLY_SUM_GENBANK
						: AssemblySummaryReader.ASSEMLY_SUM_REFSEQ));
	}
	
	@Override
	public boolean isAllowTransitiveClean() {
		return false;
	}

	@Override
	public List<File> getFiles() {
		return files;
	}

	@Override
	protected String getFTPDir(File file) {
		return getProject().getConfig().isUseGenBank() ? FTP_DIR_GENBANK : FTP_DIR_REFSEQ;
	}
}
