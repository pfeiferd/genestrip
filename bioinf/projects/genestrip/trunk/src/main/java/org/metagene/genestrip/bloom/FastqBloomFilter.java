package org.metagene.genestrip.bloom;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.nio.file.Files;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import org.metagene.genestrip.bloom.KMerBloomIndex.ByteRingBuffer;
import org.metagene.genestrip.util.CGAT;

public class FastqBloomFilter {
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

	public void runFilter(File fastgz, boolean writeIndexed, boolean writeNotIndexed) throws IOException {

		MyFileInputStream fStream = new MyFileInputStream(fastgz);
//		BufferedInputStream bStream = new BufferedInputStream(fStream, bufferSize);
		GZIPInputStream gStream = new GZIPInputStream(fStream, bufferSize);

		String name = fastgz.getName();
		int dot = name.lastIndexOf(".");
		String indexedName;
		String notIndexedName;
		if (dot != -1) {
			indexedName = name.substring(0, dot) + "_indexed" + name.substring(dot, name.length());
			notIndexedName = name.substring(0, dot) + "_not_indexed" + name.substring(dot, name.length());
		} else {
			indexedName = name + "_indexed";
			notIndexedName = name + "_not_indexed";
		}

		BufferedOutputStream bIndexed = null;
		if (writeIndexed) {
			FileOutputStream indexed = new FileOutputStream(new File(indexedName));
			GZIPOutputStream gIndexed = new GZIPOutputStream(indexed, bufferSize, false);
			bIndexed = new BufferedOutputStream(gIndexed, bufferSize);
		}
		BufferedOutputStream bNotIndexed = null;
		if (writeNotIndexed) {
			FileOutputStream notIndexed = new FileOutputStream(new File(notIndexedName));
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
			MyFileInputStream fStream) throws IOException {
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
				}
				else {
					bite = buffer[count];
				}
				lReadBuffer[line][lc[line]++] = bite;
				if (bite == '\n') {
					line++;
					if (line == 2) {
						res = classifyRead(lReadBuffer[1], lc[1] - 1);
					} else if (line == 4) {
						line = 0;
						if (res) {
							out = indexed;							
						}
						else {
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
						if (total % 1000 == 0) {
							double ratio = fStream.bRead / (double) fastqFileSize;
							long stopTime = System.currentTimeMillis();

							double diff = (stopTime - startTime);
							double totalTime = diff / ratio;
							double totalHours = totalTime / 1000 / 60 / 60;

							System.out.println("Elapse hours:" + diff / 1000 / 60 / 60);
							System.out.println("Estimated total hours:" + totalHours);
							System.out.println("Reads processed: " + total);
							System.out.println("Not indexed: " + notIndexedC);
							System.out.println("Not indexed ratio:" + ((double) notIndexedC) / total);
						}
					}
				}
			}
		}
	}

	public boolean classifyRead(byte[] read, int readSize) {
		int counter = 0;
		int negCounter = 0;
		int startAt = index.getK() - 1;
		int max = readSize - startAt;
		int posCounterThrehold = minPosCount;
		if (posCounterThrehold <= 0) {
			posCounterThrehold = (int) (max * positiveRatio);
		}
		int negCounterThreshold = max - posCounterThrehold;

		for (int i = 0; i < readSize; i++) {
			byteRingBuffer.put(read[i]);
			if (i >= startAt) {
				if (index.contains(byteRingBuffer)) {
					counter++;
					if (counter >= posCounterThrehold || counter >= minPosCount) {
						return true;
					}
				} else {
					negCounter++;
					if (negCounter > negCounterThreshold) {
						return false;
					}
				}
			}
		}

		return false;
	}

	public static class MyFileInputStream extends FileInputStream {
		private int bRead;

		public MyFileInputStream(File file) throws FileNotFoundException {
			super(file);
		}

		@Override
		public int read() throws IOException {
			bRead++;
			return super.read();
		}

		@Override
		public int read(byte[] b) throws IOException {
			bRead += b.length;
			return super.read(b);
		}

		@Override
		public int read(byte[] b, int off, int len) throws IOException {
			bRead += len;
			return super.read(b, off, len);
		}

		public int getbRead() {
			return bRead;
		}
	}
}
