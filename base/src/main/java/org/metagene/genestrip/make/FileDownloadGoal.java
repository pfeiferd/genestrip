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
import org.metagene.genestrip.io.StreamProvider;

/**
 * A {@link FileGoal} that produces its output files by downloading them over HTTP or FTP from a
 * remote base location. When an expected MD5 check sum is known for a file, it is verified after
 * download (and successful verification may be cached in a marker file), and downloads are retried
 * on check-sum failure. Subclasses supply the base URLs, the remote directory layout and the
 * download policy.
 *
 * @param <P> the type of {@link Project} this goal belongs to
 */
public abstract class FileDownloadGoal<P extends Project> extends FileGoal<P> {
	/**
	 * Configuration key controlling whether a successful check-sum verification may be recorded in an
	 * on-disk marker file to avoid recomputation.
	 */
	public static final ConfigKey CHECK_SUM_CACHE_FILE = new ConfigKey() {
		@Override
		public String getName() {
			return "checkSumCacheFile";
		}

		@Override
		public ConfigParamInfo<?> getInfo() {
			return new ConfigParamInfo.BooleanConfigParamInfo(true);
		}
	};

	/**
	 * File name suffix of the marker file used to record a successful MD5 check-sum verification.
	 */
	public static final String CHECK_SUM_OK_SUFFIX = ".md5ok";

	private final Map<File, Boolean> checkSumOkMap;
	private FTPClient ftpClient;

	/**
	 * Creates a new file download goal.
	 *
	 * @param project the project this goal belongs to
	 * @param key     the key identifying this goal
	 * @param deps    the goals this goal depends on
	 */
	@SafeVarargs
	public FileDownloadGoal(P project, GoalKey key, Goal<P>... deps) {
		super(project, key, deps);
		checkSumOkMap = new HashMap<File, Boolean>();
	}

	/**
	 * Creates the FTP client used for downloads; may be overridden.
	 *
	 * @return a new FTP client
	 */
	public FTPClient createFTPClient() {
		return new FTPClient();
	}

	/**
	 * Returns the FTP client currently used for downloads, or {@code null} if none has been created.
	 *
	 * @return the current FTP client, or {@code null}
	 */
	protected FTPClient getFtpClient() {
		return ftpClient;
	}

	/**
	 * Disconnects the given FTP client if it is still connected.
	 *
	 * @param ftpClient the FTP client to close
	 * @throws IOException if disconnecting fails
	 */
	protected void closeFTPClient(FTPClient ftpClient) throws IOException {
		if (ftpClient.isConnected()) {
			ftpClient.disconnect();
		}
	}

	/**
	 * Connects the given FTP client to the FTP base URL, logs in and switches to binary transfer
	 * mode.
	 *
	 * @param ftpClient the FTP client to connect and configure
	 * @throws IOException if connecting, logging in or configuring fails
	 */
	protected void connectLoginAndConfigure(FTPClient ftpClient) throws IOException {
		if (getLogger().isDebugEnabled()) {
			getLogger().debug("FTP Connect " + getFTPBaseURL());
		}
		ftpClient.connect(getFTPBaseURL());
		login(ftpClient);
		ftpClient.setFileType(FTP.BINARY_FILE_TYPE);
	}

	/**
	 * Logs the given FTP client in using {@link #getLogin()} and {@link #getPassword()}.
	 *
	 * @param ftpClient the FTP client to log in
	 * @throws IOException if the login fails
	 */
	protected void login(FTPClient ftpClient) throws IOException {
		if (getLogger().isDebugEnabled()) {
			getLogger().debug("FTP Login " + getLogin());
		}
		ftpClient.login(getLogin(), getPassword());
	}

	/**
	 * Returns the FTP login name; {@code "anonymous"} by default.
	 *
	 * @return the FTP login name
	 */
	protected String getLogin() {
		return "anonymous";
	}

	/**
	 * Returns the FTP password; {@code "anonymous"} by default.
	 *
	 * @return the FTP password
	 */
	protected String getPassword() {
		return "anonymous";
	}

	/**
	 * Downloads the given file from the FTP server, connecting and logging in first if not yet
	 * connected.
	 *
	 * @param ftpClient the FTP client to download with
	 * @param file      the file to download and save to
	 * @throws IOException           if the download fails
	 * @throws IllegalStateException if the file could not be fully downloaded
	 */
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

	/**
	 * Downloads the given file over HTTP, verifying its MD5 check sum on the fly when an expected
	 * check sum is available.
	 *
	 * @param file the file to download and save to
	 * @throws IOException if the download fails
	 */
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

	/**
	 * Computes and returns the hex-encoded MD5 check sum of the given file.
	 *
	 * @param file the file to compute the check sum for
	 * @return the hex-encoded MD5 check sum
	 */
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

	/**
	 * Returns the expected MD5 check sum for the given file, or {@code null} if none is known and no
	 * verification should take place; {@code null} by default.
	 *
	 * @param file the file to return the expected check sum for
	 * @return the expected hex-encoded MD5 check sum, or {@code null}
	 */
	protected String getMD5CheckSum(File file) {
		return null;
	}

	/**
	 * Builds the HTTP URL to download the given file from.
	 *
	 * @param file the file to build the URL for
	 * @return the HTTP download URL
	 */
	protected String buildHttpURL(File file) {
		return getHttpBaseURL() + getFTPDir(file) + "/" + file.getName();
	}

	/**
	 * Builds the FTP URL to download the given file from.
	 *
	 * @param file the file to build the URL for
	 * @return the FTP download URL
	 */
	protected String buildFTPURL(File file) {
		return getFTPBaseURL() + getFTPDir(file) + "/" + file.getName();
	}

	/**
	 * Returns the base URL used for HTTP downloads.
	 *
	 * @return the HTTP base URL
	 */
	protected abstract String getHttpBaseURL();

	/**
	 * Returns the host / base URL used for FTP downloads.
	 *
	 * @return the FTP base URL
	 */
	protected abstract String getFTPBaseURL();

	/**
	 * Returns whether HTTP (rather than FTP) is used for downloads.
	 *
	 * @return {@code true} if HTTP is used, {@code false} for FTP
	 */
	protected abstract boolean isUseHttp();

	/**
	 * Returns whether missing or unverifiable files are tolerated (and merely logged) instead of
	 * causing the make to fail.
	 *
	 * @return {@code true} if missing files are ignored
	 */
	protected abstract boolean isIgnoreMissingFiles();

	/**
	 * Returns the maximum number of download attempts per file before giving up.
	 *
	 * @return the maximum number of download attempts
	 */
	protected abstract int getMaxDownloadTries();

	@Override
	protected void startMake() {
		super.startMake();
		if (!isUseHttp()) {
			ftpClient = createFTPClient();
		}
	}

	/**
	 * Returns whether the given file exists as an ordinary made file and, if an expected check sum is
	 * known, matches it.
	 */
	@Override
	public boolean isMade(File file) {
		if (!super.isMade(file)) {
			return false;
		}
		return isCheckSumOk(file, null);
	}

	/**
	 * Returns whether the given file matches its expected MD5 check sum. A file without an expected
	 * check sum is always OK. If {@code computedCheckSum} is {@code null}, a cached result, an
	 * on-disk marker file or a freshly computed check sum is used as needed.
	 *
	 * @param file             the file to verify
	 * @param computedCheckSum the already computed check sum, or {@code null} to determine it as needed
	 * @return {@code true} if the file matches its expected check sum (or has none)
	 */
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

	/**
	 * Returns whether an on-disk marker file may be used to remember a successful check-sum
	 * verification, avoiding recomputation.
	 *
	 * @return {@code true} if a check-sum marker file may be used
	 */
	protected boolean isAllowCheckSumOkCacheFile() {
		return booleanConfigValue(CHECK_SUM_CACHE_FILE);
	}

	/**
	 * Returns the marker file used to record a successful check-sum verification for the given file.
	 *
	 * @param file the file to return the marker file for
	 * @return the check-sum marker file
	 */
	protected File checkSumOkCacheFile(File file) {
		return new File(file.getAbsolutePath() + CHECK_SUM_OK_SUFFIX);
	}

	/**
	 * Deletes the given file together with its check-sum marker file.
	 */
	@Override
	protected void deleteFile(File file) {
		super.deleteFile(file);
		checkSumOkCacheFile(file).delete();
	}

	/**
	 * Downloads the given file, retrying up to {@link #getMaxDownloadTries()} times while the check
	 * sum is wrong. A missing file is tolerated (and logged) when {@link #isIgnoreMissingFiles()} is
	 * {@code true}, otherwise it causes a {@link FileNotFoundException}.
	 */
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

	/**
	 * Returns whether the given file is an "additional" file that is fetched via
	 * {@link #additionalDownload(File)} rather than through the standard HTTP/FTP path; {@code false}
	 * by default.
	 *
	 * @param file the file to check
	 * @return {@code true} if the file is an additional file
	 */
	protected boolean isAdditionalFile(File file) {
		return false;
	}

	/**
	 * Downloads an additional file not covered by the standard HTTP/FTP scheme; must be overridden to
	 * be usable.
	 *
	 * @param file the additional file to download
	 * @throws FileNotFoundException always, unless overridden
	 */
	public void additionalDownload(File file) throws IOException {
		throw new FileNotFoundException("No implementation for additional download of file: " + file);
	}

	/**
	 * Returns the remote directory (relative to the base URL) that contains the given file.
	 *
	 * @param file the file to return the remote directory for
	 * @return the remote directory relative to the base URL
	 */
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
