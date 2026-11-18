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
import java.net.HttpURLConnection;
import java.net.URL;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.HashMap;
import java.util.Map;
import java.util.zip.GZIPInputStream;

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

	/** File name suffix of the file an HTTP download is written to until it is complete. */
	public static final String PART_SUFFIX = ".part";

	/** How long to wait for a connection, in milliseconds. */
	private static final int CONNECT_TIMEOUT_MS = 60 * 1000;
	/**
	 * How long a connection may deliver nothing before it counts as lost and the download is resumed,
	 * in milliseconds. Without it a stalled connection would hang the build for good.
	 */
	private static final int READ_TIMEOUT_MS = 5 * 60 * 1000;

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
	 * <p>
	 * The data go to {@link #partFile(File)} first and are moved to {@code file} only once they are
	 * complete, so that a download cut short - by a lost connection or by a killed process - never
	 * leaves a file under the final name, which {@link #isMade(File)} would take for done where no
	 * check sum is known (the release catalog has none). A connection can end early without an
	 * error, so completeness is checked, and a download that falls short is resumed from where it
	 * stopped by an HTTP range request rather than started again, up to {@link #getMaxDownloadTries()}
	 * times. {@code If-Range} makes sure the remainder comes from the same version of the remote file;
	 * a server that sends the whole file instead is read from the start again. A part left by an
	 * earlier process is not resumed, since nothing says which version of the remote file it belongs
	 * to.
	 * <p>
	 * Complete means as long as the server announced, in {@code Content-Length} or, since every request
	 * is a range request, in {@code Content-Range}. Where it announces neither, a gzip file is complete
	 * when its stream decompresses to the end ({@link #isCompleteGzip(File)}, which costs one pass over
	 * it), and any other file when the connection ended without an error - which is all there is to go
	 * by then.
	 *
	 * @param file the file to download and save to
	 * @throws IOException if the download fails, or is still incomplete after the last try
	 */
	protected void httpDownload(File file) throws IOException {
		String url = buildHttpURL(file);
		if (getLogger().isDebugEnabled()) {
			getLogger().debug("HTTP download " + url);
			getLogger().debug("Saving file " + file.toString());
		}
		File part = partFile(file);
		part.delete();
		MessageDigest messageDigest = null;
		if (getMD5CheckSum(file) != null) {
			try {
				messageDigest = MessageDigest.getInstance("MD5");
			} catch (NoSuchAlgorithmException e) {
				throw new RuntimeException(e);
			}
		}
		// The length the server announced for the whole file, -1 where it announced none, and what
		// identifies the version of the file the bytes so far belong to.
		long total = -1;
		String validator = null;
		int max = getMaxDownloadTries();
		for (int i = 1; ; i++) {
			long have = part.length();
			HttpURLConnection connection = (HttpURLConnection) new URL(url).openConnection();
			connection.setConnectTimeout(CONNECT_TIMEOUT_MS);
			connection.setReadTimeout(READ_TIMEOUT_MS);
			// A range from the start too: the answer to one names the whole length in Content-Range,
			// which some servers state where they send no Content-Length.
			connection.setRequestProperty("Range", "bytes=" + have + "-");
			if (have > 0 && validator != null) {
				connection.setRequestProperty("If-Range", validator);
			}
			try {
				int code = connection.getResponseCode();
				if (code == HttpURLConnection.HTTP_NOT_FOUND) {
					throw new FileNotFoundException(url);
				}
				if (code == 416) {
					// Range not satisfiable: nothing lies at or after the byte asked for.
					if (have == total) {
						// The last byte is in already.
						break;
					}
					if (have == 0 && total < 0) {
						// Not even a first byte: the remote file is empty.
						part.createNewFile();
						break;
					}
					// The part is longer than the remote file, so it does not belong to it.
					part.delete();
					throw new IOException("Range not satisfiable at " + have + " for " + url);
				}
				String contentRange = connection.getHeaderField("Content-Range");
				if (code == HttpURLConnection.HTTP_PARTIAL) {
					if (rangeStart(contentRange) != have) {
						part.delete();
						throw new IOException("Unexpected Content-Range " + contentRange + " for " + url);
					}
					long announced = rangeTotal(contentRange);
					if (announced >= 0) {
						total = announced;
					}
				} else if (code == HttpURLConnection.HTTP_OK) {
					total = connection.getContentLengthLong();
				} else {
					throw new IOException("HTTP status " + code + " for " + url);
				}
				boolean append = code == HttpURLConnection.HTTP_PARTIAL && have > 0;
				if (!append) {
					validator = connection.getHeaderField("ETag");
					if (validator == null) {
						validator = connection.getHeaderField("Last-Modified");
					}
					if (messageDigest != null) {
						messageDigest.reset();
					}
				}
				try (InputStream in = connection.getInputStream();
						OutputStream out = new FileOutputStream(part, append)) {
					byte[] buffer = new byte[1 << 16];
					int n;
					while ((n = in.read(buffer)) != -1) {
						out.write(buffer, 0, n);
						if (messageDigest != null) {
							messageDigest.update(buffer, 0, n);
						}
					}
				}
				long got = part.length();
				if (total < 0) {
					if (!isGzip(file) || isCompleteGzip(part)) {
						break;
					}
					throw new IOException("Received " + got + " bytes of an incomplete gzip stream");
				}
				if (got == total) {
					break;
				}
				if (got > total) {
					// More than the whole file: the parts do not fit together. Start over.
					part.delete();
				}
				throw new IOException("Received " + got + " of " + total + " bytes");
			} catch (FileNotFoundException e) {
				part.delete();
				throw e;
			} catch (IOException e) {
				if (i >= max) {
					part.delete();
					throw new IOException("Download of " + url + " incomplete after " + max + " tries", e);
				}
				if (getLogger().isWarnEnabled()) {
					getLogger().warn("Try " + i + " of " + max + ": download of " + url + " stopped at "
							+ part.length() + (total < 0 ? "" : " of " + total) + " bytes (" + e.getMessage()
							+ "), resuming.");
				}
			} finally {
				connection.disconnect();
			}
		}
		Files.move(part.toPath(), file.toPath(), StandardCopyOption.REPLACE_EXISTING);
		if (messageDigest != null) {
			// Gotta do it here since result will be cached:
			boolean res = isCheckSumOk(file, Hex.encodeHexString(messageDigest.digest()));
			if (getLogger().isDebugEnabled()) {
				getLogger().debug((res ? "OK" : "Wrong") + " MD5 check sum for file: " + file);
			}
		}
	}

	/**
	 * Returns the first byte position of a {@code Content-Range} header such as
	 * {@code bytes 100-199/1000}, or -1 where there is none.
	 *
	 * @param contentRange the header value, or {@code null}
	 * @return the first byte position, or -1
	 */
	static long rangeStart(String contentRange) {
		if (contentRange == null || !contentRange.startsWith("bytes ")) {
			return -1;
		}
		int dash = contentRange.indexOf('-');
		try {
			return dash < 0 ? -1 : Long.parseLong(contentRange.substring(6, dash).trim());
		} catch (NumberFormatException e) {
			return -1;
		}
	}

	/**
	 * Returns the whole length of a {@code Content-Range} header such as {@code bytes 100-199/1000},
	 * or -1 where it is not stated ({@code *}) or there is no header.
	 *
	 * @param contentRange the header value, or {@code null}
	 * @return the whole length, or -1
	 */
	static long rangeTotal(String contentRange) {
		if (contentRange == null) {
			return -1;
		}
		int slash = contentRange.lastIndexOf('/');
		try {
			return slash < 0 ? -1 : Long.parseLong(contentRange.substring(slash + 1).trim());
		} catch (NumberFormatException e) {
			return -1;
		}
	}

	/**
	 * Returns whether the given file is gzip-compressed, judged by its name.
	 *
	 * @param file the file
	 * @return whether its name ends with {@code .gz}
	 */
	protected boolean isGzip(File file) {
		return file.getName().endsWith(".gz");
	}

	/**
	 * Returns whether the given file holds a gzip stream that decompresses to its end. A stream cut
	 * short ends without its trailer, which shows as an exception; any other defect counts as
	 * incomplete too. A cut exactly between two members of a multi-member file goes unnoticed.
	 *
	 * @param file the file to check
	 * @return whether the gzip stream is complete
	 */
	protected boolean isCompleteGzip(File file) {
		try (InputStream in = new GZIPInputStream(Files.newInputStream(file.toPath()), 1 << 16)) {
			byte[] buffer = new byte[1 << 16];
			while (in.read(buffer) != -1) {
			}
			return true;
		} catch (IOException e) {
			return false;
		}
	}

	/**
	 * Returns the file an HTTP download of the given file is written to until it is complete.
	 *
	 * @param file the file to download
	 * @return the file holding the download in progress
	 */
	protected File partFile(File file) {
		return new File(file.getAbsolutePath() + PART_SUFFIX);
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
