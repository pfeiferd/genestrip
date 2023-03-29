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
import java.nio.file.Files;
import java.util.HashMap;
import java.util.Map;
import java.util.zip.GZIPInputStream;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.metagene.genestrip.util.ByteArrayToString;
import org.metagene.genestrip.util.ByteCountingFileInputStream;
import org.metagene.genestrip.util.CGAT;
import org.metagene.genestrip.util.CountingDigitTrie;

public class FastqTrieClassifier {
	protected final Log logger = LogFactory.getLog(getClass());

	private final KMerTrie<String> trie;
	private final byte[][] readBuffer;
	private final int[] c;
	private final int bufferSize;

	public FastqTrieClassifier(KMerTrie<String> trie, int maxReadSize) {
		bufferSize = 4096 * 8;
		this.trie = trie;
		readBuffer = new byte[4][maxReadSize];
		c = new int[4];
	}

	public Map<String, Long> runClassifier(File fastgz) throws IOException {
		if (logger.isInfoEnabled()) {
			logger.info("Reading file " + fastgz);
		}
		
		ByteCountingFileInputStream fStream = new ByteCountingFileInputStream(fastgz);
		GZIPInputStream gStream = new GZIPInputStream(fStream, bufferSize);

		long fastqFileSize = Files.size(fastgz.toPath());

		Map<String, Long> res = runClassifier(gStream, fastqFileSize, fStream);

		gStream.close();
		if (logger.isInfoEnabled()) {
			logger.info("Read file " + fastgz);
		}

		return res;
	}

	private Map<String, Long> runClassifier(InputStream fastqStream, long fastqFileSize, ByteCountingFileInputStream fStream)
			throws IOException {
		int line = 0;
		long total = 0;

		long startTime = System.currentTimeMillis();

		byte[] buffer = new byte[bufferSize];

		int size, count;
		byte bite;
		byte[][] lReadBuffer = readBuffer;
		int[] lc = c;
		long indexedC = 0;

		CountingDigitTrie root = new CountingDigitTrie();

		for (size = fastqStream.read(buffer); size != -1; size = fastqStream.read(buffer)) {
			for (count = 0; count < size; count++) {
				if (line == 2) {
					bite = CGAT.cgatToUpperCase(buffer[count]);
				} else {
					bite = buffer[count];
				}
				lReadBuffer[line][lc[line]++] = bite;
				if (bite == '\n') {
					line++;
					if (line == 2) {
						if (classifyRead(lReadBuffer[0], lReadBuffer[1], lc[1] - 1, root)) {
							indexedC++;
						}
					} else if (line == 4) {
						line = 0;
						lc[3] = lc[2] = lc[1] = lc[0] = 0;
						total++;
//						if (logger.isInfoEnabled()) {
//							if (total % 100000 == 0) {
//								double ratio = fStream.getBytesRead() / (double) fastqFileSize;
//								long stopTime = System.currentTimeMillis();
//
//								double diff = (stopTime - startTime);
//								double totalTime = diff / ratio;
//								double totalHours = totalTime / 1000 / 60 / 60;
//
//								logger.info("Elapse hours:" + diff / 1000 / 60 / 60);
//								logger.info("Estimated total hours:" + totalHours);
//								logger.info("Reads processed: " + total);
//								logger.info("Indexed: " + indexedC);
//								logger.info("Indexed ratio:" + ((double) indexedC) / total);
//							}
//						}
					}
				}
			}
		}

		Map<String, Long> counts = new HashMap<String, Long>();
		root.collect(counts);
		return counts;
	}

	public boolean classifyRead(byte[] descriptor, byte[] read, int readSize, CountingDigitTrie countingDigitTrie) {
		int max = readSize - trie.getLen() + 1;
		boolean found = false;
		for (int i = 0; i < max; i++) {
			String res = trie.get(read, i);
			if (res != null) {
				found = true;
				countingDigitTrie.inc(res);
				System.out.println(ByteArrayToString.toString(descriptor));
			}
		}
		return found;
	}
}
