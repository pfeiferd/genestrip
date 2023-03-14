package org.metagene.genestrip.gen.goals;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.zip.ZipEntry;
import java.util.zip.ZipInputStream;

import org.metagene.genestrip.gen.Config;
import org.metagene.genestrip.make.FileDownloadGoal;
import org.metagene.genestrip.tax.TaxTree;

public class TaxIdFileDownloadGoal extends FileDownloadGoal {
	public static final String TAX_DMP_ZIP = "taxdmp.zip";
	public static final String FTP_DIR = "/pub/taxonomy";

	private final File storeDir;
	private final List<File> files;

	public TaxIdFileDownloadGoal(String name, Config config) {
		super(name, config.getFtpBaseURL(), config.getHttpBaseURL(), config.isUseHttp());
		this.storeDir = config.getCommonBaseDir();
		files = new ArrayList<File>();
		files.add(new File(storeDir, TaxTree.NODES_DMP));
		files.add(new File(storeDir, TaxTree.NAMES_DMP));
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
		File zipFile = new File(storeDir, TAX_DMP_ZIP);
		if (zipFile.exists()) {
			zipFile.delete();
		}
	}

	@Override
	protected void makeFile(File file) throws IOException {
		File zipFile = new File(storeDir, TAX_DMP_ZIP);

		if (zipFile.exists() && zipFile.length() == 0) {
			if (getLogger().isInfoEnabled()) {
				getLogger().info("File delete " + zipFile.toString());
			}		
			zipFile.delete();
		}
		if (!zipFile.exists()) {
			if (isUseHttp()) {
				httpDownload(zipFile, getFTPDir(zipFile), TAX_DMP_ZIP);
			} else {
				ftpDownload(getFtpClient(), zipFile, getFTPDir(zipFile), TAX_DMP_ZIP);
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
				BufferedOutputStream out = new BufferedOutputStream(new FileOutputStream(file));
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
