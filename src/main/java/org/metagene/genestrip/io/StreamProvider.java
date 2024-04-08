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

public class StreamProvider {
	public static final int DEFAULT_BUFFER_SIZE = 8 * 4069;

	private static int bufferSize = DEFAULT_BUFFER_SIZE;

	public static void setBufferSize(int bufferSize) {
		StreamProvider.bufferSize = bufferSize;
	}

	public static int getBufferSize() {
		return bufferSize;
	}
	
	public static ByteCountingInputStreamAccess getByteCountingInputStreamForFile(File file, boolean noGZ)
			throws IOException {
		final ByteCountingFileInputStream in = new ByteCountingFileInputStream(file);
		InputStream[] res = new InputStream[1];
		res[0] = getInputStreamForFile(file, noGZ);
		if (!noGZ && isGZIPFile(file)) {
			res[0] = new GZIPInputStream(in, bufferSize);
		} else {
			res[0] = new BufferedInputStream(in, bufferSize);
		}
		return new ByteCountingInputStreamAccess() {
			@Override
			public long getBytesRead() {
				return in.getBytesRead();
			}

			@Override
			public InputStream getInputStream() {
				return res[0];
			}
		};
	}

	public static InputStream getInputStreamForFile(File file) throws IOException {
		return getInputStreamForFile(file, false);
	}

	public static InputStream getInputStreamForFile(File file, boolean noGZ) throws IOException {
		FileInputStream in = new FileInputStream(file);
		if (!noGZ && isGZIPFile(file)) {
			return new GZIPInputStream(in, bufferSize);
		} else {
			return new BufferedInputStream(in, bufferSize);
		}
	}

	public static OutputStream getOutputStreamForFile(File file) throws IOException {
		return getOutputStreamForFile(file, false);
	}

	public static OutputStream getOutputStreamForFile(File file, boolean noGZ) throws IOException {
		FileOutputStream in = new FileOutputStream(file);
		if (!noGZ && isGZIPFile(file)) {
			return new GZIPOutputStream(in, bufferSize);
		} else {
			return new BufferedOutputStream(in, bufferSize);
		}
	}

	public static boolean isGZIPFile(File file) {
		String name = file.getName();

		return name.endsWith(".gz") || name.endsWith(".gzip");
	}

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

	public static long getUncompressedSize(File gzipFile) {
		try {
			InputStream stream = StreamProvider.getInputStreamForFile(gzipFile);
			long uSize;
			for (uSize = 0; stream.read() != -1; uSize++)
				;
			return uSize;
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}

	public interface ByteCountingInputStreamAccess {
		public InputStream getInputStream();

		public long getBytesRead();
	}
}
