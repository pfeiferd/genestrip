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

import me.tongfei.progressbar.ProgressBar;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;

public class ByteCountingFileInputStream extends FileInputStream {
	private final ProgressBar progressBar;
	private long readToMark;
	private long bRead;

	public ByteCountingFileInputStream(File file) throws FileNotFoundException {
		this(file, null);
	}

	public ByteCountingFileInputStream(File file, ProgressBar progressBar) throws FileNotFoundException {
		super(file);
		bRead = 0;
		readToMark = 0;
        this.progressBar = progressBar;
    }

	@Override
	public int read() throws IOException {
		int res = super.read();
		if (res != -1) {
			bRead++;
			if (progressBar != null) {
				progressBar.stepTo(bRead);
			}
		}
		return res;
	}

	@Override
	public int read(byte[] b) throws IOException {
		int res = super.read(b);
		bRead += res;
		if (progressBar != null) {
			progressBar.stepTo(bRead);
		}
		return res;
	}

	@Override
	public int read(byte[] b, int off, int len) throws IOException {
		int res = super.read(b, off, len);
		bRead += res;
		if (progressBar != null) {
			progressBar.stepTo(bRead);
		}
		return res;
	}
	
	@Override
	public synchronized void mark(int readlimit) {
		super.mark(readlimit);
		readToMark = bRead;
	}

	@Override
	public synchronized void reset() throws IOException {
		super.reset();
		bRead = readToMark;
		if (progressBar != null) {
			progressBar.stepTo(bRead);
		}
	}
	
	public long getBytesRead() {
		return bRead;
	}
}
