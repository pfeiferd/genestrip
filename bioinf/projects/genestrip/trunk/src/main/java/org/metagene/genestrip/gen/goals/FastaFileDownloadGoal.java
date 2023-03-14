package org.metagene.genestrip.gen.goals;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.metagene.genestrip.gen.Config;
import org.metagene.genestrip.make.FileDownloadGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.tax.AssemblySummaryReader.FTPEntryWithQuality;

public class FastaFileDownloadGoal extends FileDownloadGoal {
	private final ObjectGoal<List<FTPEntryWithQuality>> entryGoal;
	private final int baseURLLen;
	private final File storeDir;

	private List<File> files;
	private Map<String, String> fileToDir;

	public FastaFileDownloadGoal(Config config, File storeDir, ObjectGoal<List<FTPEntryWithQuality>> entryGoal,
			Goal... deps) {
		super("fastasdownload", config.getFtpBaseURL(), config.getHttpBaseURL(), config.isUseHttp(), deps);
		this.entryGoal = entryGoal;
		this.storeDir = storeDir;
		baseURLLen = config.getHttpBaseURL().length();
	}

	@Override
	public List<File> getFiles() {
		if (files == null) {
			files = new ArrayList<File>();
			fileToDir = new HashMap<String, String>();
			for (FTPEntryWithQuality entry : entryGoal.get()) {
				files.add(new File(storeDir, entry.getFileName()));
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
