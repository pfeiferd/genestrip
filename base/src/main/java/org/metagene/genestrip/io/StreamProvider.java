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
package org.metagene.genestrip.io;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.Collection;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

/**
 * Utility for obtaining buffered input and output streams for files, transparently applying GZIP
 * (de)compression based on the file name, plus helpers for estimating uncompressed sizes of GZIP
 * files.
 */
public class StreamProvider {
	// Static utility class - not meant to be instantiated.
	private StreamProvider() {
	}

	/** The default buffer size used for the created streams. */
	public static final int DEFAULT_BUFFER_SIZE = 8 * 4096;

	private static int bufferSize = DEFAULT_BUFFER_SIZE;

	/**
	 * Sets the buffer size used for subsequently created streams.
	 *
	 * @param bufferSize the buffer size in bytes
	 */
	public static void setBufferSize(int bufferSize) {
		StreamProvider.bufferSize = bufferSize;
	}

	/**
	 * Returns the buffer size used for created streams.
	 *
	 * @return the buffer size in bytes
	 */
	public static int getBufferSize() {
		return bufferSize;
	}

	/**
	 * Opens a buffered input stream for the file, transparently decompressing it if its name
	 * indicates a GZIP file.
	 *
	 * @param file the file to open
	 * @return a buffered (and possibly decompressing) input stream for the file
	 * @throws IOException if the file cannot be opened or read
	 */
	public static InputStream getInputStreamForFile(File file) throws IOException {
		return getInputStreamForFile(file, false);
	}

	/**
	 * Opens a buffered input stream for the file. If {@code noGZ} is {@code true}, no GZIP
	 * decompression is applied even for {@code .gz} files.
	 *
	 * @param file the file to open
	 * @param noGZ if {@code true}, suppress GZIP decompression even for GZIP files
	 * @return a buffered (and possibly decompressing) input stream for the file
	 * @throws IOException if the file cannot be opened or read
	 */
	public static InputStream getInputStreamForFile(File file, boolean noGZ) throws IOException {
		FileInputStream in = new FileInputStream(file);
		if (!noGZ && isGZIPFile(file)) {
			return new GZIPInputStream(in, bufferSize);
		} else {
			return new BufferedInputStream(in, bufferSize);
		}
	}

	/**
	 * Opens a buffered output stream for the file, transparently GZIP-compressing what is written if
	 * its name indicates a GZIP file.
	 *
	 * @param file the file to open
	 * @return a buffered (and possibly compressing) output stream for the file
	 * @throws IOException if the file cannot be opened or written
	 */
	public static OutputStream getOutputStreamForFile(File file) throws IOException {
		return getOutputStreamForFile(file, false);
	}

	/**
	 * Opens a buffered output stream for the file. If {@code noGZ} is {@code true}, no GZIP
	 * compression is applied even for {@code .gz} files.
	 *
	 * @param file the file to open
	 * @param noGZ if {@code true}, suppress GZIP compression even for GZIP files
	 * @return a buffered (and possibly compressing) output stream for the file
	 * @throws IOException if the file cannot be opened or written
	 */
	public static OutputStream getOutputStreamForFile(File file, boolean noGZ) throws IOException {
		FileOutputStream in = new FileOutputStream(file);
		if (!noGZ && isGZIPFile(file)) {
			return new GZIPOutputStream(in, bufferSize);
		} else {
			return new BufferedOutputStream(in, bufferSize);
		}
	}

	/**
	 * Returns whether the file's name indicates a GZIP file.
	 *
	 * @param file the file whose name is inspected
	 * @return whether the file's name indicates a GZIP file
	 */
	public static boolean isGZIPFile(File file) {
		return isGZIPFileName(file.getName());
	}

	/**
	 * Returns whether the given file name ends with {@code .gz} or {@code .gzip}.
	 *
	 * @param name the file name to inspect
	 * @return whether the file name ends with {@code .gz} or {@code .gzip}
	 */
	public static boolean isGZIPFileName(String name) {
		return name.endsWith(".gz") || name.endsWith(".gzip");
	}

	/**
	 * Estimates the total uncompressed size of the given GZIP files. Up to {@code maxCheckForGuess}
	 * of the existing files are fully decompressed to measure their compression ratio, which is then
	 * extrapolated to the summed compressed size of all files. If fewer than {@code maxCheckForGuess}
	 * files exist, their measured uncompressed sizes are returned directly.
	 *
	 * @param gzipFiles the GZIP files whose combined uncompressed size is estimated
	 * @param maxCheckForGuess the maximum number of existing files to fully decompress for the estimate
	 * @return the estimated total uncompressed size of the given GZIP files
	 */
	public static long guessUncompressedSize(Collection<File> gzipFiles, int maxCheckForGuess) {
		int i = 0;
		long cSizeSum = 0;
		long uSizeSum = 0;
		for (File file : gzipFiles) {
			if (i == maxCheckForGuess) {
				break;
			}
			if (file.exists()) {
				cSizeSum += file.length();
				uSizeSum += getUncompressedSize(file);
				i++;
			}
		}
		if (i < maxCheckForGuess) {
			return uSizeSum;
		} else {
			double cRate = ((double) uSizeSum) / cSizeSum;

			long sizeSum = 0;
			for (File file : gzipFiles) {
				sizeSum += file.length();
			}
			return (long) (sizeSum * cRate);
		}
	}

	/**
	 * Returns the uncompressed size of the given GZIP file by decompressing and reading it fully.
	 *
	 * @param gzipFile the GZIP file to measure
	 * @return the uncompressed size of the given GZIP file in bytes
	 */
	public static long getUncompressedSize(File gzipFile) {
		try (InputStream stream = StreamProvider.getInputStreamForFile(gzipFile)) {
			long uSize;
			for (uSize = 0; stream.read() != -1; uSize++)
				;
			return uSize;
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}
}
