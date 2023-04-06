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
import java.io.InputStream;
import java.io.OutputStream;
import java.nio.file.Files;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.metagene.genestrip.util.CGAT;
import org.metagene.genestrip.util.CGATRingBuffer;
import org.metagene.genestrip.util.StreamProvider;
import org.metagene.genestrip.util.StreamProvider.ByteCountingInputStreamAccess;

public class FastqBloomFilter {
	protected final Log logger = LogFactory.getLog(getClass());

	private final KMerBloomIndex index;
	private final double positiveRatio;
	private final int minPosCount;
	private final CGATRingBuffer byteRingBuffer;
	private final byte[][] readBuffer;
	private final int[] c;
	private final int bufferSize;

	public FastqBloomFilter(KMerBloomIndex index, double positiveRatio, int minPosCount, int maxReadSize) {
		bufferSize = 4096 * 8;
		this.index = index;
		this.positiveRatio = positiveRatio;
		this.minPosCount = minPosCount;
		byteRingBuffer = new CGATRingBuffer(index.getK());
		readBuffer = new byte[4][maxReadSize];
		c = new int[4];
	}

	public void runFilter(File fastgz, File filteredFile, File restFile) throws IOException {
		ByteCountingInputStreamAccess access = StreamProvider.getByteCountingInputStreamForFile(restFile, false);

		OutputStream bIndexed = null;
		if (filteredFile != null) {
			bIndexed = StreamProvider.getOutputStreamForFile(filteredFile);
		}
		OutputStream bNotIndexed = null;
		if (restFile != null) {
			bNotIndexed = StreamProvider.getOutputStreamForFile(restFile);
		}

		long fastqFileSize = Files.size(fastgz.toPath());

		runFilter(access.getInputStream(), bIndexed, bNotIndexed, fastqFileSize, access);

		access.getInputStream().close();
		if (bIndexed != null) {
			bIndexed.close();
		}
		if (bNotIndexed != null) {
			bNotIndexed.close();
		}
	}

	private void runFilter(InputStream fastqStream, OutputStream indexed, OutputStream notIndexed, long fastqFileSize,
			ByteCountingInputStreamAccess byteCountAccess) throws IOException {
		int line = 0;
		boolean res = false;
		long total = 0;

		long startTime = System.currentTimeMillis();

		byte[] buffer = new byte[bufferSize];

		int size, count;
		byte bite;
		OutputStream out;
		byte[][] lReadBuffer = readBuffer;
		int[] lc = c;
		long notIndexedC = 0;

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
						res = isAcceptRead(lReadBuffer[1], lc[1] - 1);
					} else if (line == 4) {
						line = 0;
						if (res) {
							out = indexed;
						} else {
							out = notIndexed;
							notIndexedC++;
						}
						if (out != null) {
							out.write(lReadBuffer[0], 0, lc[0]);
							out.write(lReadBuffer[1], 0, lc[1]);
							out.write(lReadBuffer[2], 0, lc[2]);
							out.write(lReadBuffer[3], 0, lc[3]);
						}
						lc[3] = lc[2] = lc[1] = lc[0] = 0;
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
								logger.info("Not indexed: " + notIndexedC);
								logger.info("Not indexed ratio:" + ((double) notIndexedC) / total);
							}
						}
					}
				}
			}
		}
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
