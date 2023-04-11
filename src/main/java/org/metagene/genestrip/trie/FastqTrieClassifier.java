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

import org.metagene.genestrip.fastq.AbstractFastqReader;
import org.metagene.genestrip.util.CountingDigitTrie;
import org.metagene.genestrip.util.StreamProvider;
import org.metagene.genestrip.util.StreamProvider.ByteCountingInputStreamAccess;

public class FastqTrieClassifier extends AbstractFastqReader {
	private final KMerTrie<String> trie;
	private CountingDigitTrie root;

	private ByteCountingInputStreamAccess byteCountAccess;
	private long fastqFileSize;
	private long startTime;
	private long total;
	private long indexedC;

	public FastqTrieClassifier(KMerTrie<String> trie, int maxReadSize) {
		super(maxReadSize);
		this.trie = trie;
	}

	public Map<String, Long> runClassifier(File fastq) throws IOException {
		byteCountAccess = StreamProvider.getByteCountingInputStreamForFile(fastq, false);
		fastqFileSize = Files.size(fastq.toPath());

		startTime = System.currentTimeMillis();
		total = 0;
		indexedC = 0;

		root = new CountingDigitTrie();
		readFastq(byteCountAccess.getInputStream());

		byteCountAccess.getInputStream().close();

		Map<String, Long> counts = new HashMap<String, Long>();
		root.collect(counts);
		return counts;
	}

	@Override
	protected void nextEntry() throws IOException {
//				System.out.println(ByteArrayUtil.toString(readDescriptor));
//				System.out.println(ByteArrayUtil.toString(read));
//
		boolean res = classifyRead(false);
		res |= classifyRead(true);
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

	protected boolean classifyRead(boolean reverse) {
		int max = readSize - trie.getLen();
		boolean found = false;
		for (int i = 0; i < max; i++) {
			String res = trie.get(read, i, reverse);
			if (res != null) {
				found = true;
				root.inc(res);
			}
		}
		return found;
	}
}
