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

import com.google.common.hash.BloomFilter;
import com.google.common.hash.Funnel;
import com.google.common.hash.PrimitiveSink;

public class KMerBloomIndex implements Serializable {
	private static final long serialVersionUID = 1L;

	private final String name;
	private final Date creationDate;

	private final int k;
	private final long expectedInsertions;
	private final double expectedFpp;

	private final BloomFilter<CGATRingBuffer> index;
	private final CGATRingBuffer byteRingBuffer;

	private long n;

	private PutListener putListener;

	public KMerBloomIndex(String name, int k, long expectedInsertions, double expectedFpp, PutListener putListener) {
		this.name = name;
		this.k = k;
		this.n = 0;
		this.creationDate = new Date();
		this.expectedInsertions = expectedInsertions;
		this.expectedFpp = expectedFpp;
		this.putListener = putListener;

		index = BloomFilter.create(new MyFunnel(), expectedInsertions, expectedFpp);
		byteRingBuffer = new CGATRingBuffer(k);
	}

	public void putDirectKMer(byte[] kMer, int start) {
		byteRingBuffer.directPut = kMer;
		byteRingBuffer.directPutStart = start;
		index.put(byteRingBuffer);
	}

	public void put(byte bite) {
		byteRingBuffer.directPut = null;
		byteRingBuffer.put(bite);
		if (byteRingBuffer.filled && byteRingBuffer.isCGAT()) {
			if (putListener != null) {
				if (!index.mightContain(byteRingBuffer)) {
					putListener.newEntry(byteRingBuffer);
					index.put(byteRingBuffer);
					n++;
				} else {
					putListener.oldEntry(byteRingBuffer);
				}
			} else {
				index.put(byteRingBuffer);
				n++;
			}
		}
	}

	public void resetPut() {
		byteRingBuffer.reset();
	}

	public boolean contains(CGATRingBuffer byteRingBuffer) {
		return index.mightContain(byteRingBuffer);
	}

	public double getExpectedFpp() {
		return index.expectedFpp();
	}

	public long getExpectedInsertions() {
		return expectedInsertions;
	}

	public String getName() {
		return name;
	}

	public Date getCreationDate() {
		return creationDate;
	}

	public int getK() {
		return k;
	}

	public long getN() {
		return n;
	}

	public long getBitSize() {
		long n = expectedInsertions;
		double p = expectedFpp;
		return (long) (-n * Math.log(p) / (Math.log(2) * Math.log(2)));
	}

	private static class MyFunnel implements Funnel<CGATRingBuffer> {
		private static final long serialVersionUID = 1L;

		@Override
		public void funnel(CGATRingBuffer from, PrimitiveSink into) {
			if (from.directPut != null) {
				int max = from.data.length + from.directPutStart;
				for (int i = from.directPutStart; i < max; i++) {
					into.putByte((byte) from.directPut[i]);
				}
			} else {
				int size = from.data.length;
				for (int i = 0; i < size; i++) {
					into.putByte((byte) from.data[(from.end + i) % size]);
				}
			}
		}
	}

	public void save(File filterFile) throws IOException {
		ObjectOutputStream oOut = new ObjectOutputStream(StreamProvider.getOutputStreamForFile(filterFile));
		oOut.writeObject(this);
		oOut.close();
	}

	public static KMerBloomIndex load(File filterFile) throws IOException, ClassNotFoundException {
		ObjectInputStream oOut = new ObjectInputStream(StreamProvider.getInputStreamForFile(filterFile));
		KMerBloomIndex res = (KMerBloomIndex) oOut.readObject();
		oOut.close();
		return res;
	}

	public interface PutListener {
		public void newEntry(CGATRingBuffer buffer);

		public void oldEntry(CGATRingBuffer buffer);
	}
}
