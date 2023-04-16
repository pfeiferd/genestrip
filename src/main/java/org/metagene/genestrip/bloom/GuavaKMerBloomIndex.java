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

import org.metagene.genestrip.util.CGATRingBuffer;

import com.google.common.hash.BloomFilter;
import com.google.common.hash.Funnel;
import com.google.common.hash.PrimitiveSink;

public class GuavaKMerBloomIndex extends AbstractKMerBloomIndex {
	private static final long serialVersionUID = 1L;

	private final long expectedInsertions;
	private final double expectedFpp;

	private final BloomFilter<CGATRingBuffer> index;

	public GuavaKMerBloomIndex(String name, int k, long expectedInsertions, double expectedFpp, PutListener putListener) {
		super(name, k, putListener);
		this.expectedInsertions = expectedInsertions;
		this.expectedFpp = expectedFpp;

		index = BloomFilter.create(new MyFunnel(), expectedInsertions, expectedFpp);
	}

	public void putDirectKMer(byte[] kMer, int start) {
		byteRingBuffer.directPut = kMer;
		byteRingBuffer.directPutStart = start;
		index.put(byteRingBuffer);
	}
	
	@Override
	protected void putInternal() {
		byteRingBuffer.directPut = null;
		index.put(byteRingBuffer);
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
	
	@Override
	public int getByteSize() {
		long n = expectedInsertions;
		double p = expectedFpp;
		return (int) (-n * Math.log(p) / (Math.log(2) * Math.log(2))) / 8;
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
}
