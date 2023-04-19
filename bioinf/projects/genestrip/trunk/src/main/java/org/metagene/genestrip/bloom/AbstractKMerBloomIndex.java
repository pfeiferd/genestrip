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
package org.metagene.genestrip.bloom;

import java.io.File;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;
import java.util.Date;

import org.metagene.genestrip.util.CGATRingBuffer;
import org.metagene.genestrip.util.StreamProvider;

public abstract class AbstractKMerBloomIndex implements Serializable {
	private static final long serialVersionUID = 1L;

	private final String name;
	private final Date creationDate;
	private final PutListener putListener;

	protected final CGATRingBuffer byteRingBuffer;

	private long n;

	public AbstractKMerBloomIndex(String name, int k, PutListener putListener) {
		this.name = name;
		this.n = 0;
		this.creationDate = new Date();
		this.putListener = putListener;
		byteRingBuffer = new CGATRingBuffer(k);
	}

	public abstract void putDirectKMer(byte[] kMer, int start);

	public void put(byte bite) {
		byteRingBuffer.directPut = null;
		byteRingBuffer.put(bite);
		if (byteRingBuffer.filled && byteRingBuffer.isCGAT()) {
			if (putListener != null) {
				if (contains(byteRingBuffer)) {
					putListener.newEntry(byteRingBuffer);
					putInternal();
					n++;
				} else {
					putListener.oldEntry(byteRingBuffer);
				}
			} else {
				putInternal();
				n++;
			}
		}
	}
	
	protected abstract void putInternal();

	public void resetPut() {
		byteRingBuffer.reset();
	}

	public abstract boolean contains(CGATRingBuffer byteRingBuffer);
	
	public abstract boolean contains(byte[] seq, int start, boolean reverse);

	public abstract long getExpectedInsertions();

	public String getName() {
		return name;
	}

	public Date getCreationDate() {
		return creationDate;
	}

	public int getK() {
		return byteRingBuffer.getSize();
	}

	public long getN() {
		return n;
	}

	public abstract int getByteSize();

	public void save(File filterFile) throws IOException {
		ObjectOutputStream oOut = new ObjectOutputStream(StreamProvider.getOutputStreamForFile(filterFile));
		oOut.writeObject(this);
		oOut.close();
	}

	public static AbstractKMerBloomIndex load(File filterFile) throws IOException, ClassNotFoundException {
		ObjectInputStream oOut = new ObjectInputStream(StreamProvider.getInputStreamForFile(filterFile));
		AbstractKMerBloomIndex res = (AbstractKMerBloomIndex) oOut.readObject();
		oOut.close();
		return res;
	}

	public interface PutListener {
		public void newEntry(CGATRingBuffer buffer);

		public void oldEntry(CGATRingBuffer buffer);
	}
}