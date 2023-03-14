package org.metagene.genestrip.bloom;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.zip.GZIPInputStream;

import org.metagene.genestrip.util.CGAT;

public class FastaIndexer {
	protected InputStream getFastaGZStream(File file) throws IOException {
		FileInputStream fStream = new FileInputStream(file);
		GZIPInputStream gStream = new GZIPInputStream(fStream, 4096);

		return gStream;
	}

	public void fillIndex(KMerBloomIndex index, File file) throws IOException {
		InputStream inputStream = getFastaGZStream(file);
		byte[] buffer = new byte[4096 * 8];

		fillIndexHelp(index, inputStream, buffer);
	}

	private void fillIndexHelp(KMerBloomIndex index, InputStream inputStream, byte[] buffer) throws IOException {
//		long count = 0;
		byte bite;

//		long startTime = System.currentTimeMillis();
		boolean infoLine = false;
		int i;

		for (int size = inputStream.read(buffer); size != -1; size = inputStream.read(buffer)) {
			for (i = 0; i < size; i++) {
				bite = buffer[i];
				if (bite == '>') {
					infoLine = true;
				} else if (bite == '\n') {
					if (infoLine) {
						index.resetPut();
					}
					infoLine = false;
				} else if (!infoLine) {
					index.put(CGAT.CGAT_TO_UPPER_CASE[bite]);
//					count++;
//					if (count % 10000000 == 0) {
//						double ratio = (double) count / index.getExpectedInsertions();
//						long stopTime = System.currentTimeMillis();
//
//						double diff = (stopTime - startTime);
//						double totalTime = diff / ratio;
//						double totalHours = totalTime / 1000 / 60 / 60;
//
//						System.out.println("Elapse hours:" + diff / 1000 / 60 / 60);
//						System.out.println("Estimated total hours:" + totalHours);
//					}
				}
			}
		}

		inputStream.close();
	}
}