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

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.PrintStream;
import java.io.Serializable;
import java.util.Date;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

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

	private final BloomFilter<ByteRingBuffer> index;
	private final ByteRingBuffer byteRingBuffer;

	private long n;
	private long m;

	private PutListener putListener;

	public KMerBloomIndex(String name, int k, long expectedInsertions, double expectedFpp, PutListener putListener) {
		this.name = name;
		this.k = k;
		this.n = this.m = 0;
		this.creationDate = new Date();
		this.expectedInsertions = expectedInsertions;
		this.expectedFpp = expectedFpp;
		this.putListener = putListener;

		index = BloomFilter.create(new MyFunnel(), expectedInsertions, expectedFpp);
		byteRingBuffer = new ByteRingBuffer(k);
	}

	public void putDirectKMer(byte[] kMer, int start) {
		byteRingBuffer.directPut = kMer;
		byteRingBuffer.directPutStart = start;
		index.put(byteRingBuffer);
	}

	public void put(byte bite) {
		byteRingBuffer.directPut = null;
		byteRingBuffer.put(bite);
		m++;
		if (m >= k) {
			if (byteRingBuffer.isCGAT()) {
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
	}

	public void resetPut() {
		m = 0;
	}

	public boolean contains(ByteRingBuffer byteRingBuffer) {
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

	private static class MyFunnel implements Funnel<ByteRingBuffer> {
		private static final long serialVersionUID = 1L;

		@Override
		public void funnel(ByteRingBuffer from, PrimitiveSink into) {
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

	public static class ByteRingBuffer implements Serializable {
		private static final long serialVersionUID = 1L;

		private int end;
		private byte[] data;

		private byte[] directPut;
		private int directPutStart;

		private int invalidCPos;

		public ByteRingBuffer(int size) {
			data = new byte[size];
			end = 0;
			invalidCPos = 0;
		}

		public void put(byte c) {
			data[end] = c;
			if (c != 'C' && c != 'G' && c != 'A' && c != 'T') {
				invalidCPos = data.length;
			} else if (invalidCPos > 0) {
				invalidCPos--;
			}
			end = (end + 1) % data.length;
		}

		public int getSize() {
			return data.length;
		}

		public byte get(int index) {
			return data[(end + index) % data.length];
		}

		public String toString() {
			StringBuilder builder = new StringBuilder();

			for (int i = 0; i < data.length; i++) {
				builder.append((char) get(i));
			}

			return builder.toString();
		}

		public boolean isCGAT() {
			return invalidCPos == 0;
		}

		public void toStream(PrintStream stream) {
			for (int i = 0; i < data.length; i++) {
				stream.append((char) get(i));
			}
		}
	}

	public void save(File filterFile) throws IOException {
		ObjectOutputStream oOut = new ObjectOutputStream(
				new BufferedOutputStream(new GZIPOutputStream(new FileOutputStream(filterFile))));
		oOut.writeObject(this);
		oOut.close();
	}

	public static KMerBloomIndex load(File filterFile) throws IOException, ClassNotFoundException {
		ObjectInputStream oOut = new ObjectInputStream(
				new BufferedInputStream(new GZIPInputStream(new FileInputStream(filterFile))));
		KMerBloomIndex res = (KMerBloomIndex) oOut.readObject();
		oOut.close();
		return res;
	}

	public interface PutListener {
		public void newEntry(ByteRingBuffer buffer);

		public void oldEntry(ByteRingBuffer buffer);
	}
}
