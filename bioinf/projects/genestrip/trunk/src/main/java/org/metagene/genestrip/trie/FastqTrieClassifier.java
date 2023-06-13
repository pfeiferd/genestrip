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
import java.io.OutputStream;
import java.io.PrintStream;
import java.nio.file.Files;
import java.util.HashMap;
import java.util.Map;

import org.metagene.genestrip.fastq.AbstractFastqReader;
import org.metagene.genestrip.util.ByteArrayUtil;
import org.metagene.genestrip.util.CountingDigitTrie;
import org.metagene.genestrip.util.StreamProvider;
import org.metagene.genestrip.util.StreamProvider.ByteCountingInputStreamAccess;

public class FastqTrieClassifier extends AbstractFastqReader {
	private final KMerTrie<String> trie;
	private CountingDigitTrie root;

	private ByteCountingInputStreamAccess byteCountAccess;
	private long fastqFileSize;
	private long startTime;
	private long indexedC;

	private OutputStream indexed;

	private PrintStream out;

	public FastqTrieClassifier(KMerTrie<String> trie, int maxReadSize) {
		super(maxReadSize);
		this.trie = trie;
	}

	public Map<String, Long> runClassifier(File fastq, File filteredFile) throws IOException {

		byteCountAccess = StreamProvider.getByteCountingInputStreamForFile(fastq, false);
		fastqFileSize = Files.size(fastq.toPath());

		indexed = filteredFile != null ? StreamProvider.getOutputStreamForFile(filteredFile) : null;
		if (indexed != null) {
			OutputStream outs = StreamProvider
					.getOutputStreamForFile(new File(filteredFile.getParent(), filteredFile.getName() + ".out"));
			out = new PrintStream(outs);
		}

		startTime = System.currentTimeMillis();
		indexedC = 0;

		root = new CountingDigitTrie();
		readFastq(byteCountAccess.getInputStream());

		if (indexed != null) {
			indexed.close();
			out.close();
		}
		byteCountAccess.getInputStream().close();

		Map<String, Long> counts = new HashMap<String, Long>();
		root.collect(counts);
		return counts;
	}

	@Override
	protected void nextEntry() throws IOException {
		boolean res = classifyRead(false);
		res |= classifyRead(true);
		if (res) {
			indexedC++;
			if (indexed != null) {
				rewriteInput(indexed);
			}
		}
		if (logger.isInfoEnabled()) {
			if (reads % 100000 == 0) {
				double ratio = byteCountAccess.getBytesRead() / (double) fastqFileSize;
				long stopTime = System.currentTimeMillis();

				double diff = (stopTime - startTime);
				double totalTime = diff / ratio;
				double totalHours = totalTime / 1000 / 60 / 60;

				logger.info("Elapsed hours:" + diff / 1000 / 60 / 60);
				logger.info("Estimated total hours:" + totalHours);
				logger.info("Reads processed: " + reads);
				logger.info("Indexed: " + indexedC);
				logger.info("Indexed ratio:" + ((double) indexedC) / reads);
			}
		}
	}

	@Override
	protected void done() throws IOException {
		if (logger.isInfoEnabled()) {
			logger.info("Total indexed: " + indexedC);
		}
	}

	protected boolean classifyRead(boolean reverse) {
		String oldRes = null, res = null;
		int max = readSize - trie.getLen() + 1;
		boolean found = false;
		int c = 0;
		for (int i = 0; i < max; i++) {
			oldRes = res;
			res = trie.get(read, i, reverse);
			if (res != null) {
				if (out != null && !found) {
					out.print((reverse ? "Reverse " : "") + ByteArrayUtil.toString(readDescriptor));
				}
				found = true;
				root.inc(res);
				if (out != null && oldRes != res && oldRes != null) {
					out.print(" " + oldRes + ":" + c + ":" + (i + 1));
					c = 0;
				}
				c++;
			} else {
				if (out != null && c > 0) {
					out.print(" " + oldRes + ":" + c + ":" + (i + 1));
				}
				c = 0;
			}
		}
		if (out != null && c > 0) {
			out.print(" " + oldRes + ":" + c);
		}
		if (out != null && found) {
			out.println();
		}
		return found;
	}
}
