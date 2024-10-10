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

public class FastqBloomFilter extends AbstractLoggingFastqStreamer {
	private final double positiveRatio;
	private final int minPosCount;
	private final MurmurCGATBloomFilter filter;

	private OutputStream indexed;
	private OutputStream notIndexed;

	public FastqBloomFilter(MurmurCGATBloomFilter filter, int minPosCount, double positiveRatio, int initialReadSize,
			int maxQueueSize, ExecutionContext bundle, boolean withProbs) {
		super(filter.getK(), initialReadSize, maxQueueSize, bundle, withProbs);
		this.filter = filter;
		this.minPosCount = minPosCount;
		this.positiveRatio = positiveRatio;
	}

	@Override
	protected void nextEntry(ReadEntry readStruct, int index) throws IOException {
		boolean res = false;
		if (minPosCount > 0) {
			res = isAcceptReadByAbs(readStruct, true) || isAcceptReadByAbs(readStruct, false);
		} else {
			res = isAcceptReadByRatio(readStruct, true) || isAcceptReadByRatio(readStruct, false);
		}
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

	protected boolean isAcceptReadByAbs(ReadEntry readStruct, boolean reverse) {
		int counter = 0;
		int negCounter = 0;
		int max = readStruct.readSize - k + 1;
		int negThreshold = max - minPosCount;

		MyReadEntry re = (MyReadEntry) readStruct;

		for (int i = 0; i < max; i++) {
			if (filter.contains(readStruct.read, i, re.badPos, reverse)) {
				counter++;
				if (counter >= minPosCount) {
					return true;
				}
			} else {
				if (re.badPos[0] >= 0) {
					i += re.badPos[0];
				}
				negCounter++;
				if (negCounter > negThreshold) {
					return false;
				}
			}
		}

		return false;
	}

	protected boolean isAcceptReadByRatio(ReadEntry readStruct, boolean reverse) {
		int counter = 0;
		int negCounter = 0;
		int max = readStruct.readSize - k + 1;
		int posCounterThrehold = (int) (max * positiveRatio);
		int negCounterThreshold = max - posCounterThrehold;

		MyReadEntry re = (MyReadEntry) readStruct;

		for (int i = 0; i < max; i++) {
			if (filter.contains(readStruct.read, i, re.badPos, reverse)) {
				counter++;
				if (counter >= posCounterThrehold) {
					return true;
				}
			} else {
				if (re.badPos[0] >= 0) {
					i += re.badPos[0];
				}
				negCounter++;
				if (negCounter > negCounterThreshold) {
					return false;
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
