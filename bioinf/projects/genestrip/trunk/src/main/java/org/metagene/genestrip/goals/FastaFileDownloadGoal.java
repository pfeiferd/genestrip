package org.metagene.genestrip.goals;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.make.FileDownloadGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.tax.AssemblySummaryReader.FTPEntryWithQuality;

public class FastaFileDownloadGoal extends FileDownloadGoal<GSProject> {
	private final ObjectGoal<List<FTPEntryWithQuality>, GSProject> entryGoal;
	private final int baseURLLen;

	private List<File> files;
	private Map<String, String> fileToDir;

	@SafeVarargs
	public FastaFileDownloadGoal(GSProject project,
			ObjectGoal<List<FTPEntryWithQuality>, GSProject> entryGoal, Goal<GSProject>... deps) {
		super(project, "fastasdownload", project.getConfig().getFtpBaseURL(), project.getConfig().getHttpBaseURL(),
				project.getConfig().isUseHttp(), deps);
		this.entryGoal = entryGoal;
		baseURLLen = project.getConfig().getHttpBaseURL().length();
	}

	@Override
	public List<File> getFiles() {
		if (files == null) {
			files = new ArrayList<File>();
			fileToDir = new HashMap<String, String>();
			for (FTPEntryWithQuality entry : entryGoal.get()) {
				files.add(new File(getProject().getFastasDir(), entry.getFileName()));
				fileToDir.put(entry.getFileName(), getFtpDirFromURL(entry.getFtpURL()));
			}
		}
		return files;
	}

	protected String getFtpDirFromURL(String url) {
		return url.substring(baseURLLen);

	}

	@Override
	protected String getFTPDir(File file) {
		return fileToDir.get(file.getName());
	}

}
