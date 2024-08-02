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

import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;

public class ByteCountingInputStream extends InputStream {
	protected long bRead;
	private long readToMark;
	private final InputStream delegate;

	public ByteCountingInputStream(InputStream delegate) throws FileNotFoundException {
		this.delegate = delegate;
	}

	@Override
	public int read() throws IOException {
		bRead++;
		return delegate.read();
	}

	@Override
	public int read(byte[] b) throws IOException {
		bRead += b.length;
		return delegate.read(b);
	}

	@Override
	public int read(byte[] b, int off, int len) throws IOException {
		bRead += len;
		return delegate.read(b, off, len);
	}

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
		return delegate.skip(n);
	}
}
