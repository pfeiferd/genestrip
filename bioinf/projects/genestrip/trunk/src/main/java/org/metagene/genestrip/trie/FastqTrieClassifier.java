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
import java.nio.file.Files;
import java.util.HashMap;
import java.util.Map;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.metagene.genestrip.fastq.AbstractFastqReader;
import org.metagene.genestrip.util.CountingDigitTrie;
import org.metagene.genestrip.util.StreamProvider;
import org.metagene.genestrip.util.StreamProvider.ByteCountingInputStreamAccess;

public class FastqTrieClassifier {
	protected final Log logger = LogFactory.getLog(getClass());

	private final KMerTrie<String> trie;
	private final AbstractFastqReader fastqReader;
	
	private CountingDigitTrie root;
	private ByteCountingInputStreamAccess byteCountAccess;
	private long fastqFileSize;
	private long startTime;
	private long total;
	private long indexedC;

	public FastqTrieClassifier(KMerTrie<String> trie, int maxReadSize) {
		this.trie = trie;

		fastqReader = new AbstractFastqReader(maxReadSize) {
			@Override
			protected void nextEntry() throws IOException {
				boolean res = classifyRead(read, readSize - 1, false);
				res |= classifyRead(read, readSize - 1, true);						
				if (res) {
					indexedC++;
				}
				total++;
				if (logger.isInfoEnabled()) {
					if (total % 100000 == 0) {
						double ratio = byteCountAccess.getBytesRead() / (double) fastqFileSize;
						long stopTime = System.currentTimeMillis();

						double diff = (stopTime - startTime);
						double totalTime = diff / ratio;
						double totalHours = totalTime / 1000 / 60 / 60;

						logger.info("Elapsed hours:" + diff / 1000 / 60 / 60);
						logger.info("Estimated total hours:" + totalHours);
						logger.info("Reads processed: " + total);
						logger.info("Indexed: " + indexedC);
						logger.info("Indexed ratio:" + ((double) indexedC) / total);
					}
				}
			}
		};
	}

	public Map<String, Long> runClassifier(File fastq) throws IOException {
		root = new CountingDigitTrie();
		byteCountAccess = StreamProvider.getByteCountingInputStreamForFile(fastq, false);
		fastqFileSize = Files.size(fastq.toPath());
		
		startTime = System.currentTimeMillis();
		total = 0;
		indexedC = 0;

		fastqReader.readFastq(byteCountAccess.getInputStream());

		byteCountAccess.getInputStream().close();

		Map<String, Long> counts = new HashMap<String, Long>();
		root.collect(counts);
		return counts;
	}

	protected boolean classifyRead(byte[] read, int readSize, boolean reverse) {
		int max = readSize - trie.getLen() + 1;
		boolean found = false;
		for (int i = 0; i < max; i++) {
			String res = trie.get(read, i, reverse);
			if (res != null) {
				found = true;
				root.inc(res);
//				System.out.println(ByteArrayToString.toString(descriptor));
			}
		}
		return found;
	}
}
