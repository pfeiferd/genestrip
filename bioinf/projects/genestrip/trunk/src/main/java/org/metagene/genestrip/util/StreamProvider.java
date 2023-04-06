package org.metagene.genestrip.util;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
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

	public static ByteCountingInputStreamAccess getByteCountingInputStreamForFile(File file, boolean noGZ) throws IOException {
		final ByteCountingFileInputStream in = new ByteCountingFileInputStream(file);
		InputStream[] res = new InputStream[1];
		if (!noGZ &&  isGZIPFile(file)) {
			res[0] = new GZIPInputStream(in, bufferSize);
		}
		else {
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

	public interface ByteCountingInputStreamAccess {
		public InputStream getInputStream();

		public long getBytesRead();
	}
}
