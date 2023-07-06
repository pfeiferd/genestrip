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
package org.metagene.genestrip.match;

import java.io.Serializable;

import org.metagene.genestrip.fasta.AbstractFastaReader;
import org.metagene.genestrip.store.KMerStore;
import org.metagene.genestrip.util.CGATRingBuffer;

public abstract class FastaTrieCleaner<T extends Serializable> extends AbstractFastaReader {
	private final KMerStore<T> trie;
	private final CGATRingBuffer ringBuffer;

	public FastaTrieCleaner(KMerStore<T> trie, int bufferSize) {
		super(bufferSize);
		this.trie = trie;
		this.ringBuffer = new CGATRingBuffer(trie.getK());
	}

	@Override
	protected void dataLine() {
		for (int i = 0; i < size - 1; i++) {
			ringBuffer.put(target[i]);
			if (ringBuffer.filled && ringBuffer.isCGAT()) {
				T res = trie.get(ringBuffer, false);
//				System.out.println(res + ": Found " + ringBuffer);
				if (isMatchingValue(res)) {
//					System.out.println("Removing " + ringBuffer);
					trie.put(ringBuffer, null, false);
				} else {
					res = trie.get(ringBuffer, true);
//					System.out.println(res + ": Found reverse " + ringBuffer);
					if (isMatchingValue(res)) {
//						System.out.println("Removing " + ringBuffer);
						trie.put(ringBuffer, null, true);
					}
				}
			}
		}
	}

	protected abstract boolean isMatchingValue(T value);
}