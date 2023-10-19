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
package org.metagene.genestrip.make;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.net.URL;
import java.nio.channels.Channels;
import java.nio.channels.ReadableByteChannel;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.net.ftp.FTP;
import org.apache.commons.net.ftp.FTPClient;
import org.metagene.genestrip.io.StreamProvider;
import org.metagene.genestrip.make.FileDownloadGoal.DownloadProject;

public abstract class FileDownloadGoal<P extends DownloadProject> extends FileGoal<P> {
	private FTPClient ftpClient;
	private List<File> availableFiles;

	@SafeVarargs
	public FileDownloadGoal(P project, String name, Goal<P>... deps) {
		super(project, name, deps);
		availableFiles = new ArrayList<File>();
	}

	public List<File> getAvailableFiles() {
		return availableFiles;
	}

	public FTPClient createFTPClient() {
		return new FTPClient();
	}

	protected FTPClient getFtpClient() {
		return ftpClient;
	}

	protected void closeFTPClient(FTPClient ftpClient) throws IOException {
		if (ftpClient.isConnected()) {
			ftpClient.disconnect();
		}
	}

	protected void connectLoginAndConfigure(FTPClient ftpClient) throws IOException {
		if (getLogger().isInfoEnabled()) {
			getLogger().info("FTP Connect " + getProject().getBaseFTPURL());
		}
		ftpClient.connect(getProject().getBaseFTPURL());
		login(ftpClient);
		ftpClient.setFileType(FTP.BINARY_FILE_TYPE);
	}

	protected void login(FTPClient ftpClient) throws IOException {
		if (getLogger().isInfoEnabled()) {
			getLogger().info("FTP Login " + getLogin());
		}
		ftpClient.login(getLogin(), getPassword());
	}

	protected String getLogin() {
		return "anonymous";
	}

	protected String getPassword() {
		return "anonymous";
	}

	protected void ftpDownload(FTPClient ftpClient, File file) throws IOException {
		if (!ftpClient.isConnected()) {
			connectLoginAndConfigure(ftpClient);
		}

		if (getLogger().isInfoEnabled()) {
			getLogger().info("Saving file " + file.toString());
		}
		OutputStream out = StreamProvider.getOutputStreamForFile(file, true);
		ftpClient.changeWorkingDirectory(getFTPDir(file));
		if (getLogger().isInfoEnabled()) {
			getLogger().info("FTP download " + buildFTPURL(file));
		}
		if (!ftpClient.retrieveFile(file.getName(), out)) {
			out.close();
			file.delete();
			throw new IllegalStateException("Could not fully download " + file.getName());
		} else {
			out.close();
		}
	}

	protected void httpDownload(File file) throws IOException {
		String url = buildHttpURL(file);
		if (getLogger().isInfoEnabled()) {
			getLogger().info("HTTP download " + url);
		}
		ReadableByteChannel readableByteChannel = Channels.newChannel(new URL(url).openStream());
		if (getLogger().isInfoEnabled()) {
			getLogger().info("Saving file " + file.toString());
		}
		FileOutputStream out = new FileOutputStream(file);
		out.getChannel().transferFrom(readableByteChannel, 0, Long.MAX_VALUE);
		out.close();
	}

	protected String buildHttpURL(File file) {
		return getHttpBaseURL() + getFTPDir(file) + "/" + file.getName();
	}

	protected String buildFTPURL(File file) {
		return getFTPBaseURL() + getFTPDir(file) + "/" + file.getName();
	}
	
	protected String getHttpBaseURL() {
		return getProject().getHttpBaseURL();
	}

	protected String getFTPBaseURL() {
		return getProject().getBaseFTPURL();
	}
	
	@Override
	protected void startMake() {
		if (!getProject().isUseHttp()) {
			ftpClient = createFTPClient();
		}
		availableFiles.clear();
	}

	@Override
	public boolean isMade(File file) {
		boolean res = super.isMade(file);
		if (res && !availableFiles.contains(file)) {
			availableFiles.add(file);
		}
		return res;
	}

	@Override
	protected void makeFile(File file) throws IOException {
		try {
			if (isAdditionalFile(file)) {
				additionalDownload(file);
			} else if (isUseHttp()) {
				httpDownload(file);
			} else {
				ftpDownload(ftpClient, file);
			}
			if (!availableFiles.contains(file)) {
				availableFiles.add(file);
			}
		} catch (FileNotFoundException e) {
			if (!getProject().isIgnoreMissingFiles()) {
				throw e;
			}
			if (isAdditionalFile(file)) {
				getLogger().warn("Missing file for download " + file.getCanonicalPath());
			} else if (getProject().isUseHttp()) {
				getLogger().warn("Missing file for download " + buildHttpURL(file));
			} else {
				getLogger().warn("Missing file for download " + buildFTPURL(file));
			}
		}
	}
	
	protected boolean isUseHttp() {
		return getProject().isUseHttp();
	}

	protected boolean isAdditionalFile(File file) {
		return false;
	}

	public void additionalDownload(File file) throws IOException {
		throw new FileNotFoundException("No implementation for additional download of file: " + file);
	}

	protected abstract String getFTPDir(File file);

	@Override
	protected void endMake() {
		if (ftpClient != null) {
			try {
				closeFTPClient(ftpClient);
			} catch (IOException e) {
				throw new RuntimeException(e);
			}
		}
	}

	@Override
	public String toString() {
		return "download file goal: " + getName();
	}

	public interface DownloadProject {
		public boolean isUseHttp();

		public String getBaseFTPURL();

		public String getHttpBaseURL();

		public boolean isIgnoreMissingFiles();
	}
}
