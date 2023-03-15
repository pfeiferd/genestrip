package org.metagene.genestrip.util;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;

public class ByteCountingFileInputStream extends FileInputStream {
	private long bRead;

	public ByteCountingFileInputStream(File file) throws FileNotFoundException {
		super(file);
	}

	@Override
	public int read() throws IOException {
		bRead++;
		return super.read();
	}

	@Override
	public int read(byte[] b) throws IOException {
		bRead += b.length;
		return super.read(b);
	}

	@Override
	public int read(byte[] b, int off, int len) throws IOException {
		bRead += len;
		return super.read(b, off, len);
	}

	public long getBytesRead() {
		return bRead;
	}
}
