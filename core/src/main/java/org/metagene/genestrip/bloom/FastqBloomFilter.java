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

import org.metagene.genestrip.ExecutionContext;
import org.metagene.genestrip.fastq.AbstractLoggingFastqStreamer;
import org.metagene.genestrip.io.StreamProvider;
import org.metagene.genestrip.io.StreamingResourceStream;
import org.metagene.genestrip.util.CGAT;

/**
 * Streams FASTQ reads and classifies each read by how many of its (canonical) k-mers are found in a
 * {@link KMerProbFilter}: reads with enough matching k-mers are written to one output and the rest to
 * another.
 */
public class FastqBloomFilter extends AbstractLoggingFastqStreamer {
	private final double positiveRatio;
	private final int minPosCount;
	private final KMerProbFilter filter;

	private OutputStream indexed;
	private OutputStream notIndexed;

	/**
	 * Creates a filter that accepts a read when enough of its k-mers are present in {@code filter}.
	 *
	 * @param k              the k-mer length
	 * @param filter         the k-mer probability filter to test membership against
	 * @param minPosCount    if positive, the absolute number of matching k-mers required to accept a read;
	 *                       otherwise {@code positiveRatio} of the read's k-mers must match
	 * @param positiveRatio  fraction of matching k-mers required when {@code minPosCount} is not positive
	 * @param initialReadSize the initial read buffer size in bytes
	 * @param maxQueueSize   the maximum size of the read processing queue
	 * @param bundle         the execution context driving the parallel processing
	 * @param withProbs      if {@code true}, per-base quality probabilities are read and preserved
	 */
	public FastqBloomFilter(int k, KMerProbFilter filter, int minPosCount, double positiveRatio, int initialReadSize,
							int maxQueueSize, ExecutionContext bundle, boolean withProbs) {
		super(k, initialReadSize, maxQueueSize, bundle, withProbs);
		this.filter = filter;
		this.minPosCount = minPosCount;
		this.positiveRatio = positiveRatio;
	}

	/**
	 * Filters the given FASTQ streams, writing accepted reads to {@code filteredFile} and the rest to
	 * {@code restFile}; either file may be {@code null} to discard that side.
	 *
	 * @param fastqs       the FASTQ input streams to filter
	 * @param filteredFile the file to write accepted reads to, or {@code null} to discard them
	 * @param restFile     the file to write rejected reads to, or {@code null} to discard them
	 * @throws IOException if reading the input or writing the output fails
	 */
	public void runFilter(StreamingResourceStream fastqs, File filteredFile, File restFile) throws IOException {
		try (OutputStream lindexed = filteredFile != null ? StreamProvider.getOutputStreamForFile(filteredFile) : null;
			 OutputStream lnotIndexed = restFile != null ? StreamProvider.getOutputStreamForFile(restFile) : null) {
			indexed = lindexed;
			notIndexed = lnotIndexed;
			processFastqStreams(fastqs);
		}
		indexed = null;
		notIndexed = null;
	}

	@Override
	protected void nextEntry(ReadEntry readStruct, int index) throws IOException {
		MyReadEntry re = (MyReadEntry) readStruct;

		boolean res = isAcceptRead(re);
		if (res) {
			if (indexed != null) {
				rewriteInput(readStruct, indexed);
			}
		} else {
			if (notIndexed != null) {
				rewriteInput(readStruct, notIndexed);
			}
		}
	}

	@Override
	protected ReadEntry createReadEntry(int maxReadSizeBytes, boolean withProbs, Object... config) {
		return new MyReadEntry(maxReadSizeBytes, withProbs);
	}

	/**
	 * Decides whether the given read should be accepted based on how many of its canonical k-mers are
	 * present in the filter.
	 *
	 * @param entry the read to test
	 * @return whether the given read has enough of its canonical k-mers present in the filter, scanning
	 *         only until the accept or reject threshold is decided.
	 */
	protected boolean isAcceptRead(final MyReadEntry entry) {
		int max = entry.readSize - k + 1;
		int posThreshold = (minPosCount > 0) ? minPosCount : (int) (max * positiveRatio);
		int negThreshold = max - posThreshold;

		long kmer = -1;
		long reverseKmer = -1;
		int counter = 0;
		int negCounter = 0;
		for (int i = 0; i < max; i++) {
			if (kmer == -1) {
				kmer = CGAT.kMerToLongStraight(entry.read, i, k, entry.badPos);
				if (kmer == -1) {
					i = entry.badPos[0];
				} else {
					reverseKmer = CGAT.kMerToLongReverse(entry.read, i, k, null);
				}
			} else {
				kmer = CGAT.nextKMerStraight(kmer, entry.read[i + k - 1], k);
				if (kmer == -1) {
					i += k - 1;
				} else {
					reverseKmer = CGAT.nextKMerReverse(reverseKmer, entry.read[i + k - 1], k);
				}
			}
			if (kmer != -1) {
				if (filter.containsLong(CGAT.standardKMer(kmer, reverseKmer))) {
					counter++;
					if (counter >= posThreshold) {
						return true;
					}
				} else {
					negCounter++;
					if (negCounter > negThreshold) {
						return false;
					}
				}
			}
		}

		return false;
	}


	/**
	 * A read entry that additionally tracks the position of the last invalid base encountered.
	 */
	protected static class MyReadEntry extends ReadEntry {
		/** Single-element holder for the position of the last invalid base. */
		public final int[] badPos = new int[1];

		/**
		 * Creates the read entry.
		 *
		 * @param maxReadSizeBytes the maximum read size in bytes
		 * @param withProbs        if {@code true}, per-base quality probabilities are stored
		 */
		protected MyReadEntry(int maxReadSizeBytes, boolean withProbs) {
			super(maxReadSizeBytes, withProbs);
		}
	}
}
