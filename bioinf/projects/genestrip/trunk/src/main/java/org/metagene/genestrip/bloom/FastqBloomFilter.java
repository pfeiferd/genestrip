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

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.nio.file.Files;

import org.metagene.genestrip.fastq.AbstractFastqReader;
import org.metagene.genestrip.util.StreamProvider;
import org.metagene.genestrip.util.StreamProvider.ByteCountingInputStreamAccess;

public class FastqBloomFilter extends AbstractFastqReader {
	private static final byte[] LINE_3 = new byte[] { '+', '\n' };

	private final AbstractKMerBloomIndex index;
	private final int k;

	private final double positiveRatio;
	private final int minPosCount;

	private OutputStream indexed;
	private OutputStream notIndexed;
	private OutputStream out;
	private ByteCountingInputStreamAccess byteCountAccess;
	private long fastqFileSize;
	private long startTime;
	private long indexedC;

	public FastqBloomFilter(AbstractKMerBloomIndex index, int minPosCount, double positiveRatio, int maxReadSize) {
		super(maxReadSize);
		this.index = index;
		this.k = index.getK();
		this.minPosCount = minPosCount;
		this.positiveRatio = positiveRatio;
	}

	@Override
	protected void nextEntry() throws IOException {
		boolean res = false;
		if (minPosCount > 0) {
			res = isAcceptReadByAbs(true) || isAcceptReadByAbs(false);
		} else {
			res = isAcceptReadByRatio(true) || isAcceptReadByRatio(false);
		}
		if (res) {
			out = indexed;
			indexedC++;
		} else {
			out = notIndexed;
		}
		if (out != null) {
			rewriteInput(out);
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

	public void runFilter(File fastq, File filteredFile, File restFile) throws IOException {
		byteCountAccess = StreamProvider.getByteCountingInputStreamForFile(fastq, false);
		fastqFileSize = Files.size(fastq.toPath());

		indexed = filteredFile != null ? StreamProvider.getOutputStreamForFile(filteredFile) : null;
		notIndexed = restFile != null ? StreamProvider.getOutputStreamForFile(restFile) : null;

		startTime = System.currentTimeMillis();
		indexedC = 0;

		readFastq(byteCountAccess.getInputStream());

		if (indexed != null) {
			indexed.close();
		}
		if (notIndexed != null) {
			notIndexed.close();
		}

		byteCountAccess.getInputStream().close();
	}
	
	private final int[] badPos = new int[1];

	protected boolean isAcceptReadByAbs(boolean reverse) {
		int counter = 0;
		int negCounter = 0;
		int max = readSize - k + 1;
		int negThreshold = max - minPosCount;

		for (int i = 0; i < max; i++) {
			badPos[0] = -1;
			if (index.contains(read, i, reverse, badPos)) {
				counter++;
				if (counter >= minPosCount) {
					return true;
				}
			} else {
				if (badPos[0] > 0) {
					i += badPos[0];
				}
				negCounter++;
				if (negCounter > negThreshold) {
					return false;
				}
			}
		}

		return false;
	}

	protected boolean isAcceptReadByRatio(boolean reverse) {
		int counter = 0;
		int negCounter = 0;
		int max = readSize - k + 1;
		int posCounterThrehold = (int) (max * positiveRatio);
		int negCounterThreshold = max - posCounterThrehold;
		
		for (int i = 0; i < max; i++) {
			badPos[0] = -1;
			if (index.contains(read, i, reverse, badPos)) {
				counter++;
				if (counter >= posCounterThrehold) {
					return true;
				}
			} else {
				if (badPos[0] > 0) {
					i += badPos[0];
				}
				negCounter++;
				if (negCounter > negCounterThreshold) {
					return false;
				}
			}
		}

		return false;
	}

}
