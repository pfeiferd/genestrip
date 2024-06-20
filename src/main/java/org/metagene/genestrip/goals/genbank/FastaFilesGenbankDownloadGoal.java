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
import java.io.FileOutputStream;
import java.io.IOException;
import java.net.URL;
import java.nio.channels.Channels;
import java.nio.channels.ReadableByteChannel;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.metagene.genestrip.GSGoalKey;
import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.genbank.AssemblySummaryReader.FTPEntryWithQuality;
import org.metagene.genestrip.goals.GSFileDownloadGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;

public class FastaFilesGenbankDownloadGoal extends GSFileDownloadGoal {
	private final ObjectGoal<Map<TaxIdNode, List<FTPEntryWithQuality>>, GSProject> entryGoal;
	private final int baseURLLen;

	private List<File> files;
	private Map<String, Object> fileToDir;

	@SafeVarargs
	public FastaFilesGenbankDownloadGoal(GSProject project,
			ObjectGoal<Map<TaxIdNode, List<FTPEntryWithQuality>>, GSProject> entryGoal, Goal<GSProject>... deps) {
		super(project, GSGoalKey.FASTAGSENBANKDL, append(deps, entryGoal));
		this.entryGoal = entryGoal;
		baseURLLen = getHttpBaseURL().length();
	}

	@Override
	public boolean isAllowTransitiveClean() {
		return false;
	}

	@Override
	public List<File> getFiles() {
		if (files == null) {
			files = new ArrayList<File>();
			fileToDir = new HashMap<String, Object>();
			for (List<FTPEntryWithQuality> list : entryGoal.get().values()) {
				for (FTPEntryWithQuality entry : list) {
					File file = entryToFile(entry);
					if (entry.getFtpURL() != null) {
						String dir = getFtpDirFromURL(entry.getFtpURL());
						if (dir != null) {
							files.add(file);
							fileToDir.put(file.getName(), dir);
						}
					} else {
						URL url = entry.getURL();
						if (url != null) {
							// Relative path for file url is with respect to projects base directory ...
							if ("file".equals(url.getProtocol())) {
								if (!new File(url.getPath()).isAbsolute()) {
									try {
										url = new URL("file", null,
												new File(getProject().getCommon().getBaseDir(), url.getFile()).getCanonicalPath());
									} catch (IOException e) {
										throw new RuntimeException(e);
									}
								}
							}
							files.add(file);
							fileToDir.put(file.getName(), url);
						}
					}
				}
			}
			if (getLogger().isInfoEnabled()) {
				getLogger().info("Number of Fasta files to download: " + files.size());
			}
		}
		return files;
	}

	public File entryToFile(FTPEntryWithQuality entry) {
		return new File(getFastaDir(), entry.getFileName());
	}

	protected File getFastaDir() {
		return getProject().getCommon().getGenbankDir();
	}

	protected String getFtpDirFromURL(String url) {
		if (url.length() <= baseURLLen) {
			return null;
		}
		return url.substring(baseURLLen);
	}

	@Override
	protected String getFTPDir(File file) {
		return (String) fileToDir.get(file.getName());
	}

	@Override
	protected boolean isAdditionalFile(File file) {
		return fileToDir.get(file.getName()) instanceof URL;
	}

	@Override
	public void additionalDownload(File file) throws IOException {
		URL url = (URL) fileToDir.get(file.getName());

		if (getLogger().isInfoEnabled()) {
			getLogger().info("Additional download for " + url.toExternalForm());
			getLogger().info("Saving file " + file.toString());
		}
		try (ReadableByteChannel readableByteChannel = Channels.newChannel(url.openStream());
				FileOutputStream out = new FileOutputStream(file)) {
			out.getChannel().transferFrom(readableByteChannel, 0, Long.MAX_VALUE);
		}
	}
}
