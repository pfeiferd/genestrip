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
import java.io.InputStream;
import java.io.OutputStream;
import java.net.URL;
import java.nio.channels.Channels;
import java.nio.channels.ReadableByteChannel;
import java.nio.file.Files;
import java.security.DigestInputStream;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.HashMap;
import java.util.Map;

import org.apache.commons.codec.binary.Hex;
import org.apache.commons.codec.digest.DigestUtils;
import org.apache.commons.io.FileUtils;
import org.apache.commons.net.ftp.FTP;
import org.apache.commons.net.ftp.FTPClient;
import org.metagene.genestrip.GSConfigKey;
import org.metagene.genestrip.io.StreamProvider;

public abstract class FileDownloadGoal<P extends Project> extends FileGoal<P> {
	public static final String CHECK_SUM_OK_SUFFIX = ".md5ok";

	private final Map<File, Boolean> checkSumOkMap;
	private FTPClient ftpClient;

	@SafeVarargs
	public FileDownloadGoal(P project, GoalKey key, Goal<P>... deps) {
		super(project, key, deps);
		checkSumOkMap = new HashMap<File, Boolean>();
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
		if (getLogger().isDebugEnabled()) {
			getLogger().debug("FTP Connect " + getFTPBaseURL());
		}
		ftpClient.connect(getFTPBaseURL());
		login(ftpClient);
		ftpClient.setFileType(FTP.BINARY_FILE_TYPE);
	}

	protected void login(FTPClient ftpClient) throws IOException {
		if (getLogger().isDebugEnabled()) {
			getLogger().debug("FTP Login " + getLogin());
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
		if (getLogger().isDebugEnabled()) {
			getLogger().debug("Saving file " + file.toString());
			getLogger().debug("FTP download " + buildFTPURL(file));
		}
		try (OutputStream out = StreamProvider.getOutputStreamForFile(file, true)) {
			ftpClient.changeWorkingDirectory(getFTPDir(file));
			if (!ftpClient.retrieveFile(file.getName(), out)) {
				file.delete();
				throw new IllegalStateException("Could not fully download " + file.getName());
			}
		}
	}

	protected void httpDownload(File file) throws IOException {
		String url = buildHttpURL(file);
		if (getLogger().isDebugEnabled()) {
			getLogger().debug("HTTP download " + url);
			getLogger().debug("Saving file " + file.toString());
		}
		if (getMD5CheckSum(file) != null) {
			try {
				MessageDigest messageDigest = MessageDigest.getInstance("MD5");
				try (ReadableByteChannel readableByteChannel = Channels
						.newChannel(new DigestInputStream(new URL(url).openStream(), messageDigest));
						FileOutputStream out = new FileOutputStream(file)) {
					out.getChannel().transferFrom(readableByteChannel, 0, Long.MAX_VALUE);
				}
				// Gotta do it here since result will be cached:
				boolean res = isCheckSumOk(file, Hex.encodeHexString(messageDigest.digest()));
				if (getLogger().isDebugEnabled()) {
					getLogger().debug((res ? "OK" : "Wrong") + " MD5 check sum for file: " + file);
				}
			} catch (NoSuchAlgorithmException e) {
				throw new RuntimeException(e);
			}
		} else {
			try (ReadableByteChannel readableByteChannel = Channels.newChannel(new URL(url).openStream());
					FileOutputStream out = new FileOutputStream(file)) {
				out.getChannel().transferFrom(readableByteChannel, 0, Long.MAX_VALUE);
			}
		}
	}

	protected String computeMD5CheckSum(File file) {
		String res = null;
		if (getLogger().isDebugEnabled()) {
			getLogger().debug("Computing MD5 check sum for file: " + file);
		}
		try (InputStream is = Files.newInputStream(file.toPath())) {
			res = DigestUtils.md5Hex(is);
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
		return res;
	}

	protected String getMD5CheckSum(File file) {
		return null;
	}

	protected String buildHttpURL(File file) {
		return getHttpBaseURL() + getFTPDir(file) + "/" + file.getName();
	}

	protected String buildFTPURL(File file) {
		return getFTPBaseURL() + getFTPDir(file) + "/" + file.getName();
	}

	protected abstract String getHttpBaseURL();

	protected abstract String getFTPBaseURL();

	protected abstract boolean isUseHttp();

	protected abstract boolean isIgnoreMissingFiles();

	protected abstract int getMaxDownloadTries();

	@Override
	protected void startMake() {
		super.startMake();
		if (!isUseHttp()) {
			ftpClient = createFTPClient();
		}
	}

	@Override
	public boolean isMade(File file) {
		if (!super.isMade(file)) {
			return false;
		}
		return isCheckSumOk(file, null);
	}

	protected boolean isCheckSumOk(File file, String computedCheckSum) {
		String originalCheckSum = getMD5CheckSum(file);
		if (originalCheckSum == null) {
			return true;
		}
		// The map is a cache for efficiency reasons, as check sum (re-)computation can
		// be expensive.
		if (computedCheckSum == null) {
			Boolean check = checkSumOkMap.get(file);
			if (check != null) {
				return check;
			}
			if (isAllowCheckSumOkCacheFile()) {
				if (checkSumOkCacheFile(file).exists()) {
					return true;
				}
			}
			computedCheckSum = computeMD5CheckSum(file);
		}
		boolean check = originalCheckSum.equals(computedCheckSum);
		checkSumOkMap.put(file, check);
		if (check && isAllowCheckSumOkCacheFile()) {
			try {
				FileUtils.touch(checkSumOkCacheFile(file));
			} catch (IOException e) {
				// Ignore on puropse.
			}
		}
		return check;
	}

	protected boolean isAllowCheckSumOkCacheFile() {
		return booleanConfigValue(GSConfigKey.CHECK_SUM_CACHE_FILE);
	}

	protected File checkSumOkCacheFile(File file) {
		return new File(file.getAbsolutePath() + CHECK_SUM_OK_SUFFIX);
	}

	@Override
	protected void deleteFile(File file) {
		super.deleteFile(file);
		checkSumOkCacheFile(file).delete();
	}

	@Override
	protected void makeFile(File file) throws IOException {
		boolean checkSumWrong = true;
		try {
			int max = getMaxDownloadTries();
			for (int i = 0; i < max && checkSumWrong; i++) {
				if (file.exists()) {
					file.delete();
				}
				// Invalidate checkSum chache entry before download....
				checkSumOkMap.remove(file);
				if (isAdditionalFile(file)) {
					additionalDownload(file);
				} else if (isUseHttp()) {
					httpDownload(file);
				} else {
					ftpDownload(ftpClient, file);
				}
				checkSumWrong = !isCheckSumOk(file, null);
				if (checkSumWrong) {
					getLogger().warn("Trie " + (i + 1) + " of " + max + ": Bad check sum for downloaded file " + file);
				}
			}
			if (checkSumWrong) {
				throw new FileNotFoundException("Bad check sum for downloaded file " + file);
			}
		} catch (FileNotFoundException e) {
			if (!isIgnoreMissingFiles()) {
				throw e;
			}
			else if (getLogger().isWarnEnabled()) {
				if (checkSumWrong) {
					getLogger().warn("Bad check sum for downloaded file " + file.getCanonicalPath());
				} else if (isAdditionalFile(file)) {
					getLogger().warn("Missing file for download " + file.getCanonicalPath());
				} else if (isUseHttp()) {
					getLogger().warn("Missing file for download " + buildHttpURL(file));
				} else {
					getLogger().warn("Missing file for download " + buildFTPURL(file));
				}
			}
		}
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
		return "download file goal: " + getKey().getName();
	}
}
