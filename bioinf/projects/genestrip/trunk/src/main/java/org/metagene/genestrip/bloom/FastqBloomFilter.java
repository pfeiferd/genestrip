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

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.metagene.genestrip.fastq.AbstractFastqReader;
import org.metagene.genestrip.util.CGATRingBuffer;
import org.metagene.genestrip.util.StreamProvider;
import org.metagene.genestrip.util.StreamProvider.ByteCountingInputStreamAccess;

public class FastqBloomFilter {
	protected final Log logger = LogFactory.getLog(getClass());

	private static final byte[] LINE_3 = new byte[] { '+', '\n' };

	private final KMerBloomIndex index;
	private final double positiveRatio;
	private final int minPosCount;
	private final CGATRingBuffer byteRingBuffer;

	private final AbstractFastqReader fastqReader;

	private OutputStream indexed;
	private OutputStream notIndexed;
	private OutputStream out;
	private ByteCountingInputStreamAccess byteCountAccess;
	private long fastqFileSize;
	private long startTime;
	private long total;
	private long indexedC;

	public FastqBloomFilter(KMerBloomIndex index, double positiveRatio, int minPosCount, int maxReadSize) {
		this.index = index;
		this.positiveRatio = positiveRatio;
		this.minPosCount = minPosCount;
		byteRingBuffer = new CGATRingBuffer(index.getK());

		fastqReader = new AbstractFastqReader(maxReadSize) {
			@Override
			protected void nextEntry() throws IOException {
				boolean res = isAcceptRead(read, readSize - 1);
				if (res) {
					out = indexed;
					indexedC++;
				} else {
					out = notIndexed;
				}
				if (out != null) {
					out.write(readDescriptor, 0, readDescriptorSize);
					out.write(read, 0, readSize);
					out.write(LINE_3);
					out.write(readProbs, 0, readProbsSize);
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

	public void runFilter(File fastq, File filteredFile, File restFile) throws IOException {
		byteCountAccess = StreamProvider.getByteCountingInputStreamForFile(fastq, false);
		this.fastqFileSize = Files.size(fastq.toPath());

		this.indexed = filteredFile != null ? StreamProvider.getOutputStreamForFile(filteredFile) : null;
		this.notIndexed = restFile != null ? StreamProvider.getOutputStreamForFile(restFile) : null;

		startTime = System.currentTimeMillis();
		total = 0;
		indexedC = 0;

		fastqReader.readFastq(byteCountAccess.getInputStream());
		
		if (indexed != null) {
			indexed.close();
		}
		if (notIndexed != null) {
			notIndexed.close();
		}
		
		byteCountAccess.getInputStream().close();
	}

	public boolean isAcceptRead(byte[] read, int readSize) {
		int counter = 0;
		int negCounter = 0;
		int max = readSize - index.getK() - 1;
		int posCounterThrehold = 0;
		if (minPosCount <= 0) {
			posCounterThrehold = (int) (max * positiveRatio);
		}
		int negCounterThreshold = max - posCounterThrehold;

		for (int i = 0; i < readSize; i++) {
			byteRingBuffer.put(read[i]);
			if (byteRingBuffer.filled) {
				if (index.contains(byteRingBuffer)) {
					counter++;
					if (minPosCount > 0) {
						if (counter >= minPosCount) {
							return true;
						}
					} else if (counter >= posCounterThrehold) {
						return true;
					}
				} else {
					negCounter++;
					if (minPosCount > 0) {
						if (counter >= max - minPosCount) {
							return true;
						}
					} else if (negCounter > negCounterThreshold) {
						return false;
					}
				}
			}
		}

		return false;
	}
}