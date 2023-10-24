package org.metagene.genestrip.accuracy;

import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.compress.archivers.ArchiveEntry;
import org.apache.commons.compress.archivers.tar.TarArchiveInputStream;
import org.apache.commons.compress.compressors.gzip.GzipCompressorInputStream;
import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.make.FileDownloadGoal;
import org.metagene.genestrip.make.Goal;

public class AccuracyDataDownloadGoal extends FileDownloadGoal<GSProject> {
	public static final String ACCURACY_FILE = "accuracy.tgz";
	public static final String TIMING_FILE = "timing.tgz";
	public static final String FTP_DIR = "/software/kraken/dl";

	private final List<File> files;

	@SafeVarargs
	public AccuracyDataDownloadGoal(GSProject project, String name, Goal<GSProject>... deps) {
		super(project, name, deps);
		files = new ArrayList<File>();
		files.add(new File(getProject().getFastaDir(), ACCURACY_FILE));
		files.add(new File(getProject().getFastaDir(), TIMING_FILE));
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
		return FTP_DIR;
	}

	@Override
	protected void makeFile(File file) throws IOException {
		File gzipTarFile = file;

		if (gzipTarFile.exists() && gzipTarFile.length() == 0) {
			if (getLogger().isInfoEnabled()) {
				getLogger().info("File delete " + gzipTarFile.toString());
			}
			gzipTarFile.delete();
		}
		if (!gzipTarFile.exists()) {
			httpDownload(gzipTarFile);
		}

		if (getLogger().isInfoEnabled()) {
			getLogger().info("File extract " + gzipTarFile.toString());
		}
		BufferedInputStream inputStream = new BufferedInputStream(new FileInputStream(gzipTarFile));
		TarArchiveInputStream tar = new TarArchiveInputStream(new GzipCompressorInputStream(inputStream));
		ArchiveEntry entry;
		while ((entry = tar.getNextEntry()) != null) {
			Path extractTo = getProject().getFastaDir().toPath().resolve(entry.getName());
			if (entry.isDirectory()) {
				Files.createDirectories(extractTo);
			} else {
				Files.copy(tar, extractTo);
			}
		}
	}

	@Override
	protected String getHttpBaseURL() {
		return "http://ccb.jhu.edu/";
	}

	@Override
	protected String getFTPBaseURL() {
		throw new IllegalStateException("Should never be called.");
	}
}
