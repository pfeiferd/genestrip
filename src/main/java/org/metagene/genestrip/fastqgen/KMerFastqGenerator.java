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
package org.metagene.genestrip.fastqgen;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.List;
import java.util.zip.GZIPOutputStream;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.metagene.genestrip.bloom.KMerBloomIndex;
import org.metagene.genestrip.bloom.KMerBloomIndex.PutListener;
import org.metagene.genestrip.util.CGATRingBuffer;
import org.metagene.genestrip.bloom.FastaIndexer;

public class KMerFastqGenerator {
	private final Log logger = LogFactory.getLog(KMerFastqGenerator.class);

	private final int kmerSize;

	public KMerFastqGenerator(int kmerSize) {
		this.kmerSize = kmerSize;
	}

	protected Log getLogger() {
		return logger;
	}

	public int getKmerSize() {
		return kmerSize;
	}

	public long run(List<File> fastaFiles, File fastq, String name, long expectedSize) throws IOException {
		if (getLogger().isInfoEnabled()) {
			getLogger().info("File create " + fastq.toString());
		}
		FileOutputStream fStream = new FileOutputStream(fastq);
		GZIPOutputStream gStream = new GZIPOutputStream(fStream, 4096);
		FastQWriter fastQWriter = new FastQWriter(name, gStream);
		KMerBloomIndex index = new KMerBloomIndex(name, getKmerSize(), expectedSize, 0.0001, fastQWriter);
		
		if (getLogger().isInfoEnabled()) {
			getLogger().info("Bloom filter array size of " + index + ": " + index.getBitSize() / 8 / 1024 + "KB");
		}

		FastaIndexer fastaIndexer = new FastaIndexer();
		byte[] buffer = new byte[4096 * 8];

		for (File fasta : fastaFiles) {
			if (fasta.exists()) {
				if (getLogger().isInfoEnabled()) {
					getLogger().info("File read " + fasta.toString());
				}
				fastaIndexer.fillIndex(index, fasta, buffer);
			} else {
				throw new IllegalStateException("Missing fasta file " + fasta);
			}
		}

		fastQWriter.close();
		
		return fastQWriter.getAdded();
	}

	public static class FastQWriter implements PutListener {
		private final Log logger = LogFactory.getLog(FastQWriter.class);

		private final String id;
		private final PrintStream printStream;
		private long added;
		private long hits;

		public FastQWriter(String id, OutputStream out) throws IOException {
			this.id = id;
			printStream = new PrintStream(out);
			printHeader();
		}

		public long getAdded() {
			return added;
		}

		protected Log getLogger() {
			return logger;
		}

		@Override
		public void newEntry(CGATRingBuffer buffer) {
			added++;
			printBeforeRead();
			buffer.toStream(printStream);
			printStream.println();
			printAfterRead(buffer.getSize());
		}

		@Override
		public void oldEntry(CGATRingBuffer buffer) {
			hits++;
			if (getLogger().isInfoEnabled()) {
				if (hits % 100000 == 0) {
					String hitsMessage = "Hits: " + hits + " Hits Ratio: " + (double) hits / (hits + added);
					logger.info(hitsMessage);
				}
			}
		}

		public void close() {
			printStream.close();
			if (getLogger().isInfoEnabled()) {
				String hitsMessage = "Total Hits: " + hits + " Hits Ratio: " + (double) hits / (hits + added);
				logger.info(hitsMessage);
			}
		}

		protected void printHeader() {
		}

		protected void printBeforeRead() {
			printBasicID();
			printStream.print(id);
			printStream.print(':');
			printStream.println(added);
		}

		protected void printBasicID() {
			printStream.print("@GENESTRIP:");
		}

		protected void printAfterRead(int readLength) {
			printStream.println("+");
			for (int i = 0; i < readLength; i++) {
				printStream.print('~');
			}
			printStream.println();
		}
	}
}
