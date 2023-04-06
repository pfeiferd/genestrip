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
package org.metagene.genestrip.trie;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.Serializable;

import org.metagene.genestrip.util.CGATRingBuffer;
import org.metagene.genestrip.util.StreamProvider;

public abstract class FastaTrieCleaner<T extends Serializable> {
	public void cleanTrie(KMerTrie<T> trie, File file, byte[] buffer, CGATRingBuffer ringBuffer)
			throws IOException {
		InputStream inputStream = StreamProvider.getInputStreamForFile(file);
		cleanTrieHelp(trie, inputStream, buffer, ringBuffer);
		inputStream.close();
	}

	// value must not be null...
	public void cleanTrieHelp(KMerTrie<T> trie, InputStream inputStream, byte[] buffer,
			CGATRingBuffer ringBuffer) throws IOException {
		if (ringBuffer.getSize() != trie.getLen()) {
			throw new IllegalArgumentException("trie and ring buffer must have equal size");
		}
		byte bite;

		boolean infoLine = false;
		int i;

		for (int size = inputStream.read(buffer); size != -1; size = inputStream.read(buffer)) {
			for (i = 0; i < size; i++) {
				bite = buffer[i];
				if (bite == '>') {
					infoLine = true;
				} else if (bite == '\n') {
					infoLine = false;
					if (infoLine) {
						ringBuffer.reset();
					}
				} else if (!infoLine) {
					ringBuffer.put(bite);
					if (ringBuffer.filled && ringBuffer.isCGAT()) {
						T res = trie.get(ringBuffer, false);
						if (res != null) {
							System.out.println(res + ": Found " + ringBuffer);							
						}
						if (isMatchingValue(res)) {
							System.out.println("Removing " + ringBuffer);
							trie.put(ringBuffer, null, false);
						}
						else {
							res = trie.get(ringBuffer, true);
							if (res != null) {
								System.out.println(res + ": Found reverse " + ringBuffer);							
							}
							if (isMatchingValue(res)) {
								System.out.println("Removing " + ringBuffer);
								trie.put(ringBuffer, null, true);
							}
						}
					}
				}
			}
		}
	}
	
	protected abstract boolean isMatchingValue(T value);
}