package org.metagene.genestrip.make;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.net.URL;
import java.nio.channels.Channels;
import java.nio.channels.ReadableByteChannel;

import org.apache.commons.net.ftp.FTP;
import org.apache.commons.net.ftp.FTPClient;

public abstract class FileDownloadGoal extends FileGoal {
	private final String baseFTPURL;
	private final String httpBaseURL;
	private final boolean useHttp;

	private FTPClient ftpClient;

	public FileDownloadGoal(String name, String baseFTPURL, String httpBaseURL, boolean useHttp, Goal... deps) {
		super(name, deps);
		this.baseFTPURL = baseFTPURL;
		this.httpBaseURL = httpBaseURL;
		this.useHttp = useHttp;
	}

	public boolean isUseHttp() {
		return useHttp;
	}

	protected String getBaseFTPURL() {
		return baseFTPURL;
	}

	public String getHttpBaseURL() {
		return httpBaseURL;
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
			getLogger().info("FTP Connect " + getBaseFTPURL());
		}		
		ftpClient.connect(getBaseFTPURL());
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

	protected void ftpDownload(FTPClient ftpClient, File file, String ftpFolder, String fileName) throws IOException {
		if (!ftpClient.isConnected()) {
			connectLoginAndConfigure(ftpClient);
		}

		if (getLogger().isInfoEnabled()) {
			getLogger().info("File save " + file.toString());
		}		
		OutputStream out = new BufferedOutputStream(new FileOutputStream(file));
		ftpClient.changeWorkingDirectory(ftpFolder);
		if (getLogger().isInfoEnabled()) {
			getLogger().info("FTP download " + baseFTPURL + ftpFolder + "/" + fileName);
		}		
		if (!ftpClient.retrieveFile(fileName, out)) {
			out.close();
			file.delete();
			throw new IllegalStateException("Could not fully download " + fileName);
		} else {
			out.close();
		}
	}

	protected void httpDownload(File file, String path, String fileName) throws IOException {
		String url = buildHttpURL(path, fileName);
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

	protected String buildHttpURL(String path, String file) {
		return getHttpBaseURL() + path + "/" + file;
	}

	@Override
	protected void startMake() {
		if (!isUseHttp()) {
			ftpClient = createFTPClient();
		}
	}

	@Override
	protected void makeFile(File file) throws IOException {
		if (isUseHttp()) {
			httpDownload(file, getFTPDir(file), file.getName());
		} else {
			ftpDownload(ftpClient, file, getFTPDir(file), file.getName());
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
		return "download file goal: "+ getName();
	}
}
