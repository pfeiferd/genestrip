package org.metagene.genestrip.gen.goals;

import java.io.File;
import java.util.Collections;
import java.util.List;

import org.metagene.genestrip.gen.Project;
import org.metagene.genestrip.make.FileDownloadGoal;
import org.metagene.genestrip.tax.AssemblySummaryReader;

public class AssemblyFileDownloadGoal extends FileDownloadGoal<Project> {
	public static final String FTP_DIR = "/genomes/ASSEMBLY_REPORTS";
	
	private final List<File> files;
	
	public AssemblyFileDownloadGoal(Project project) {
		super(project, "assemblydownload", project.getConfig().getFtpBaseURL(), project.getConfig().getHttpBaseURL(), project.getConfig().isUseHttp());
		files = Collections.singletonList(new File(project.getConfig().getCommonBaseDir(), AssemblySummaryReader.ASSEMLY_SUM));
	}
	
	@Override
	protected List<File> getFiles() {
		return files;
	}
	
	@Override
	protected String getFTPDir(File file) {
		return FTP_DIR;
	}
}
