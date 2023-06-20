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

import java.io.EOFException;
import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.Collection;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.metagene.genestrip.bloom.AbstractKMerBloomIndex.PutListener;
import org.metagene.genestrip.bloom.FastaIndexer;
import org.metagene.genestrip.bloom.GenestripKMerBloomIndex;
import org.metagene.genestrip.util.CGATRingBuffer;
import org.metagene.genestrip.util.StreamProvider;

public class KMerFastqGenerator {
	private final Log logger = LogFactory.getLog(KMerFastqGenerator.class);

	private final FastQWriter fastQWriter;
	private final FastaIndexer fastaIndexer;
	private final GenestripKMerBloomIndex index;
	private final boolean ignoreMissingFastas;

	public KMerFastqGenerator(String name, int kmerSize, long expectedSize, double fpp, int readBufferSize,
			boolean ignoreMissingFastas) {
		fastQWriter = new FastQWriter(name);
		index = new GenestripKMerBloomIndex(name, kmerSize, expectedSize, fpp, fastQWriter);
		fastaIndexer = new FastaIndexer(index, readBufferSize);
		this.ignoreMissingFastas = ignoreMissingFastas;
	}

	public void setOutputStream(OutputStream out) {
		fastQWriter.setOutputStream(out);
	}

	protected Log getLogger() {
		return logger;
	}

	public long add(Collection<File> fastaFiles) throws IOException {
		long expectedSize = StreamProvider.guessUncompressedSize(fastaFiles, 10);
		if (getLogger().isInfoEnabled()) {
			getLogger().info("Estimated total uncompressed size: " + (expectedSize / 1024 / 1024) + " MB");
		}
		index.getFilter().clearAndEnsureCapacity(expectedSize);
		if (getLogger().isInfoEnabled()) {
			getLogger().info("Bloom filter array size of " + index + ": " + index.getFilter().getByteSize() / 1024 + " KB");
		}

		int max = fastaFiles.size();
		int counter = 0;
		for (File fasta : fastaFiles) {
			counter++;
			if (fasta.exists()) {
				if (getLogger().isInfoEnabled()) {
					getLogger().info("Processing fasta file (" + counter + "/" + max + "):" + fasta);
				}
				try {
					fastaIndexer.readFasta(fasta);
				} catch (EOFException e) {
					if (getLogger().isErrorEnabled()) {
						getLogger().error("Error when reading fasta file " + fasta, e);
					}
				}
			} else {
				if (getLogger().isErrorEnabled()) {
					getLogger().error("Missing fasta file " + fasta);
				}
				if (!ignoreMissingFastas) {
					throw new IllegalStateException("Missing fasta file " + fasta);
				}
			}
		}

		return fastQWriter.getAdded();
	}

	public static class FastQWriter implements PutListener {
		private final Log logger = LogFactory.getLog(FastQWriter.class);

		private final String id;
		private PrintStream printStream;
		private long added;
		private long hits;

		public FastQWriter(String id) {
			this.id = id;
		}

		public void setOutputStream(OutputStream out) {
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
					logger.info("Hits: " + hits + " Hits Ratio: " + (double) hits / (hits + added));
				}
			}
		}

		public void done() {
			if (getLogger().isInfoEnabled()) {
				logger.info("Total Hits: " + hits + " Hits Ratio: " + (double) hits / (hits + added));
				logger.info("Total added reads: " + added);
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
