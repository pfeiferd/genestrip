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
import java.io.IOException;
import java.nio.file.Files;
import java.text.MessageFormat;
import java.util.ArrayList;
import java.util.List;

import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.make.FileGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.util.ArraysUtil;

public class RefSeqCatalogDownloadGoal extends RefSeqDownloadGoal {
	public static final String CATALOG_FOLDER = RELEASE_FOLDER + "/release-catalog";
	public static final String CATALOG_BASE_NAME = "RefSeq-release{0}.catalog.gz";
	public static final String FILES_INSTALLED_NAME = "release{0}.files.installed";

	private final List<File> files;

	@SafeVarargs
	public RefSeqCatalogDownloadGoal(GSProject project, String name, FileGoal<GSProject> releaseNumberGoal,
			Goal<GSProject>... deps) {
		super(project, name, ArraysUtil.append(deps, releaseNumberGoal));

		try {
			byte[] encoded = Files.readAllBytes(releaseNumberGoal.getFile().toPath());
			String releaseNumber = new String(encoded).trim();

			String catalog = MessageFormat.format(CATALOG_BASE_NAME, releaseNumber);
			String filesInstalled = MessageFormat.format(FILES_INSTALLED_NAME, releaseNumber);

			files = new ArrayList<File>();

			files.add(new File(project.getConfig().getRefSeqDir(), catalog));
			files.add(new File(project.getConfig().getRefSeqDir(), filesInstalled));
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}

	public File getCatalogFile() {
		return files.get(0);
	}

	public File getInstalledFilesFile() {
		return files.get(1);
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
		return CATALOG_FOLDER;
	}
}
