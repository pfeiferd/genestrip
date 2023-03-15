package org.metagene.genestrip.bloom;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.nio.file.Files;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.metagene.genestrip.bloom.KMerBloomIndex.ByteRingBuffer;
import org.metagene.genestrip.util.ByteCountingFileInputStream;
import org.metagene.genestrip.util.CGAT;

public class FastqBloomFilter {
	protected final Log logger = LogFactory.getLog(getClass());

	private final KMerBloomIndex index;
	private final double positiveRatio;
	private final int minPosCount;
	private final ByteRingBuffer byteRingBuffer;
	private final byte[][] readBuffer;
	private final int[] c;
	private final int bufferSize;

	public FastqBloomFilter(KMerBloomIndex index, double positiveRatio, int minPosCount, int maxReadSize) {
		bufferSize = 4096 * 8;
		this.index = index;
		this.positiveRatio = positiveRatio;
		this.minPosCount = minPosCount;
		byteRingBuffer = new ByteRingBuffer(index.getK());
		readBuffer = new byte[4][maxReadSize];
		c = new int[4];
	}

	public void runFilter(File fastgz, File filteredFile, File restFile) throws IOException {
		ByteCountingFileInputStream fStream = new ByteCountingFileInputStream(fastgz);
		GZIPInputStream gStream = new GZIPInputStream(fStream, bufferSize);

		BufferedOutputStream bIndexed = null;
		if (filteredFile != null) {
			FileOutputStream indexed = new FileOutputStream(filteredFile);
			GZIPOutputStream gIndexed = new GZIPOutputStream(indexed, bufferSize, false);
			bIndexed = new BufferedOutputStream(gIndexed, bufferSize);
		}
		BufferedOutputStream bNotIndexed = null;
		if (restFile != null) {
			FileOutputStream notIndexed = new FileOutputStream(restFile);
			GZIPOutputStream gNotIndexed = new GZIPOutputStream(notIndexed, bufferSize, false);
			bNotIndexed = new BufferedOutputStream(gNotIndexed, bufferSize);
		}

		long fastqFileSize = Files.size(fastgz.toPath());

		runFilter(gStream, bIndexed, bNotIndexed, fastqFileSize, fStream);

		gStream.close();
		if (bIndexed != null) {
			bIndexed.close();
		}
		if (bNotIndexed != null) {
			bNotIndexed.close();
		}
	}

	private void runFilter(InputStream fastqStream, OutputStream indexed, OutputStream notIndexed, long fastqFileSize,
			ByteCountingFileInputStream fStream) throws IOException {
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
					bite = CGAT.CGAT_TO_UPPER_CASE[buffer[count]];
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
							if (total % 1000 == 0) {
								double ratio = fStream.getBytesRead() / (double) fastqFileSize;
								long stopTime = System.currentTimeMillis();

								double diff = (stopTime - startTime);
								double totalTime = diff / ratio;
								double totalHours = totalTime / 1000 / 60 / 60;

								logger.info("Elapse hours:" + diff / 1000 / 60 / 60);
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
		int startAt = index.getK() - 1;
		int max = readSize - startAt;
		int posCounterThrehold = 0;
		if (minPosCount <= 0) {
			posCounterThrehold = (int) (max * positiveRatio);
		}
		int negCounterThreshold = max - posCounterThrehold;

		for (int i = 0; i < readSize; i++) {
			byteRingBuffer.put(read[i]);
			if (i >= startAt) {
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
