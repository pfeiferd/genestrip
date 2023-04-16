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

import org.metagene.genestrip.bloom.KMerBloomIndex.PutListener;
import org.metagene.genestrip.util.CGATRingBuffer;

public class GenestripKMerBloomIndex extends AbstractKMerBloomIndex  {
	private static final long serialVersionUID = 1L;

	private final CGATBloomFilter index;

	public GenestripKMerBloomIndex(String name, int k, long expectedInsertions, double expectedFpp, PutListener putListener) {
		super(name, k, expectedInsertions, expectedFpp, putListener);
		index = new CGATBloomFilter(k, expectedInsertions, expectedFpp);
	}

	public void putDirectKMer(byte[] kMer, int start) {
		index.put(kMer, start);
	}

	public boolean contains(CGATRingBuffer byteRingBuffer) {
		return index.contains(byteRingBuffer, false);
	}

	public long getExpectedInsertions() {
		return index.getExpectedInsertions();
	}

	public int getK() {
		return index.getK();
	}
	
	@Override
	public int getByteSize() {
		return index.getByteSize();
	}
	
	@Override
	protected void putInternal() {
		index.put(byteRingBuffer);
	}
}