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
import java.io.FileInputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.List;
import java.util.zip.ZipEntry;
import java.util.zip.ZipInputStream;

import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.make.FileDownloadGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.tax.TaxTree;
import org.metagene.genestrip.util.StreamProvider;

public class TaxIdFileDownloadGoal extends FileDownloadGoal<GSProject> {
	public static final String TAX_DMP_ZIP = "taxdmp.zip";
	public static final String FTP_DIR = "/pub/taxonomy";

	private final List<File> files;

	@SafeVarargs
	public TaxIdFileDownloadGoal(GSProject project, String name, Goal<GSProject>... deps) {
		super(project, name, deps);
		files = new ArrayList<File>();
		files.add(new File(getProject().getConfig().getCommonBaseDir(), TaxTree.NODES_DMP));
		files.add(new File(getProject().getConfig().getCommonBaseDir(), TaxTree.NAMES_DMP));
	}

	@Override
	protected List<File> getFiles() {
		return files;
	}

	@Override
	protected String getFTPDir(File file) {
		return FTP_DIR;
	}

	@Override
	public void cleanThis() {
		super.cleanThis();
		File zipFile = new File(getProject().getConfig().getCommonBaseDir(), TAX_DMP_ZIP);
		if (zipFile.exists()) {
			zipFile.delete();
		}
	}

	@Override
	protected void makeFile(File file) throws IOException {
		File zipFile = new File(getProject().getConfig().getCommonBaseDir(), TAX_DMP_ZIP);

		if (zipFile.exists() && zipFile.length() == 0) {
			if (getLogger().isInfoEnabled()) {
				getLogger().info("File delete " + zipFile.toString());
			}
			zipFile.delete();
		}
		if (!zipFile.exists()) {
			if (getProject().isUseHttp()) {
				httpDownload(zipFile);
			} else {
				ftpDownload(getFtpClient(), zipFile);
			}
		}

		if (getLogger().isInfoEnabled()) {
			getLogger().info("File extract " + zipFile.toString());
		}
		ZipInputStream zis = new ZipInputStream(new FileInputStream(zipFile));

		byte[] buffer = new byte[1024];
		for (ZipEntry zipEntry = zis.getNextEntry(); zipEntry != null; zipEntry = zis.getNextEntry()) {
			String entryName = zipEntry.getName();
			if (entryName.equals(file.getName())) {
				if (getLogger().isInfoEnabled()) {
					getLogger().info("File save " + file.toString());
				}
				OutputStream out = StreamProvider.getOutputStreamForFile(file);
				int len;
				while ((len = zis.read(buffer)) > 0) {
					out.write(buffer, 0, len);
				}
				out.close();
			}
		}

		zis.closeEntry();
		zis.close();
	}
}
