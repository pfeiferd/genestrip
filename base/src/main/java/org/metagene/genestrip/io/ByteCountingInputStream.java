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

import java.io.IOException;
import java.io.InputStream;

/**
 * An {@link InputStream} that wraps a delegate stream and counts the number of bytes read from it,
 * keeping the count consistent across {@code mark}/{@code reset} and {@code skip}.
 */
public class ByteCountingInputStream extends InputStream {
	/** The number of bytes read from the delegate stream so far. */
	protected long bRead;
	private long readToMark;
	private final InputStream delegate;

	/**
	 * Creates a byte-counting stream wrapping the given delegate stream.
	 *
	 * @param delegate the underlying input stream to read from
	 */
	public ByteCountingInputStream(InputStream delegate) {
		this.delegate = delegate;
	}

	@Override
	public int read() throws IOException {
		int b = delegate.read();
		if (b != -1) {
			bRead++;
		}
		return b;
	}

	@Override
	public int read(byte[] b) throws IOException {
		int n = delegate.read(b);
		if (n > 0) {
			bRead += n;
		}
		return n;
	}

	@Override
	public int read(byte[] b, int off, int len) throws IOException {
		int n = delegate.read(b, off, len);
		if (n > 0) {
			bRead += n;
		}
		return n;
	}

	/**
	 * Returns the total number of bytes read so far.
	 *
	 * @return the number of bytes read
	 */
	public long getBytesRead() {
		return bRead;
	}
	
	@Override
	public void close() throws IOException {
		delegate.close();
	}

	@Override
	public int available() throws IOException {
		return delegate.available();
	}

	@Override
	public synchronized void mark(int readlimit) {
		delegate.mark(readlimit);
		readToMark = bRead;
	}

	@Override
	public boolean markSupported() {
		return delegate.markSupported();
	}

	@Override
	public synchronized void reset() throws IOException {
		delegate.reset();
		bRead = readToMark;
	}

	@Override
	public long skip(long n) throws IOException {
		long skipped = delegate.skip(n);
		if (skipped > 0) {
			bRead += skipped;
		}
		return skipped;
	}
}
