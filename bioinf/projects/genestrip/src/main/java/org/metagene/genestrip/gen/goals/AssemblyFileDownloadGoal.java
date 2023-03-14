package org.metagene.genestrip.gen.goals;

import java.io.File;
import java.util.Collections;
import java.util.List;

import org.metagene.genestrip.gen.Config;
import org.metagene.genestrip.make.FileDownloadGoal;
import org.metagene.genestrip.tax.AssemblySummaryReader;

public class AssemblyFileDownloadGoal extends FileDownloadGoal {
	public static final String FTP_DIR = "/genomes/ASSEMBLY_REPORTS";
	
	private final List<File> files;
	
	public AssemblyFileDownloadGoal(Config config) {
		super("assemblydownload", config.getFtpBaseURL(), config.getHttpBaseURL(), config.isUseHttp());
		files = Collections.singletonList(new File(config.getCommonBaseDir(), AssemblySummaryReader.ASSEMLY_SUM));
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
