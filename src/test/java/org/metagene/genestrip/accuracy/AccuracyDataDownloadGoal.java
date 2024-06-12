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
import org.metagene.genestrip.goals.GSFileDownloadGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.GoalKey.DefaultGoalKey;

public class AccuracyDataDownloadGoal extends GSFileDownloadGoal {
	public static final String ACCURACY_FILE = "accuracy.tgz";
	public static final String TIMING_FILE = "timing.tgz";
	public static final String FTP_DIR = "/software/kraken/dl";

	private final List<File> files;

	@SafeVarargs
	public AccuracyDataDownloadGoal(GSProject project, Goal<GSProject>... deps) {
		super(project, new DefaultGoalKey("accdatadownload"), deps);
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
