package org.metagene.genestrip.make;

import java.io.BufferedOutputStream;
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
			getLogger().info("File save " + file.toString());
		}
		OutputStream out = new BufferedOutputStream(new FileOutputStream(file));
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
			getLogger().info("File save " + file.toString());
		}
		FileOutputStream out = new FileOutputStream(file);
		out.getChannel().transferFrom(readableByteChannel, 0, Long.MAX_VALUE);
		out.close();
	}

	protected String buildHttpURL(File file) {
		return getProject().getHttpBaseURL() + getFTPDir(file) + "/" + file.getName();
	}

	protected String buildFTPURL(File file) {
		return getProject().getBaseFTPURL() + getFTPDir(file) + "/" + file.getName();
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
			if (getProject().isUseHttp()) {
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
			if (getProject().isUseHttp()) {
				getLogger().warn("Missing file for download " + buildHttpURL(file));
			} else {
				getLogger().warn("Missing file for download " + buildFTPURL(file));
			}
		}
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
