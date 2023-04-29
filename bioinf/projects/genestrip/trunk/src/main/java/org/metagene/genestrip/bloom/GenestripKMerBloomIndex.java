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

public class GenestripKMerBloomIndex extends AbstractKMerBloomIndex {
	private static int MAX_ARRAY_SIZE = Integer.MAX_VALUE - 8;

	private static final long serialVersionUID = 1L;

	private final AbstractCGATBloomFilter index;

	public GenestripKMerBloomIndex(String name, int k, long expectedInsertions, double expectedFpp,
			PutListener putListener) {
		super(name, k, putListener);
		index = expectedInsertions > MAX_ARRAY_SIZE ? new LargeMurmurCGATBloomFilter(k, expectedInsertions, expectedFpp)
				: new MurmurCGATBloomFilter(k, expectedInsertions, expectedFpp);
	}

	@Override
	public void putDirectKMer(byte[] kMer, int start) {
		index.put(kMer, start, null);
	}

	@Override
	public boolean contains(CGATRingBuffer byteRingBuffer, int[] badPos) {
		return index.contains(byteRingBuffer, false, badPos);
	}

	@Override
	public boolean contains(byte[] seq, int start, boolean reverse, int[] badPos) {
		return index.contains(seq, start, reverse, badPos);
	}

	@Override
	public long getExpectedInsertions() {
		return index.getExpectedInsertions();
	}

	@Override
	public long getByteSize() {
		return index.getByteSize();
	}

	@Override
	protected void putInternal() {
		index.put(byteRingBuffer, null);
	}
}