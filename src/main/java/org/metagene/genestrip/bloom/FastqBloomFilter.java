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

public class FastqBloomFilter extends AbstractLoggingFastqStreamer {
	private final double positiveRatio;
	private final int minPosCount;
	private final AbstractKMerBloomFilter filter;
	private final boolean inlined;

	private OutputStream indexed;
	private OutputStream notIndexed;

	public FastqBloomFilter(AbstractKMerBloomFilter filter, int minPosCount, double positiveRatio, int initialReadSize,
							int maxQueueSize, ExecutionContext bundle, boolean withProbs, boolean inlined) {
		super(filter.getK(), initialReadSize, maxQueueSize, bundle, withProbs);
		this.filter = filter;
		this.minPosCount = minPosCount;
		this.positiveRatio = positiveRatio;
		this.inlined = inlined;
	}

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
				}
				else {
					reverseKmer = CGAT.kMerToLongReverse(entry.read, i, k, null);
				}
			} else {
				kmer = CGAT.nextKMerStraight(kmer, entry.read[i + k - 1], k);
				if (kmer == -1) {
					i += k - 1;
				}
				else {
					reverseKmer = CGAT.nextKMerReverse(kmer, entry.read[i + k - 1], k);
				}
			}
			if (kmer != -1) {
				if (filter.containsViaHash(CGAT.standardKMer(kmer, reverseKmer))) {
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


	protected static class MyReadEntry extends ReadEntry {
		public final int[] badPos = new int[1];

		protected MyReadEntry(int maxReadSizeBytes, boolean withProbs) {
			super(maxReadSizeBytes, withProbs);
		}
	}
}
