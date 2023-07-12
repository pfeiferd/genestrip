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
package org.metagene.genestrip.kraken;

import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.metagene.genestrip.util.BufferedLineReader;
import org.metagene.genestrip.util.ByteArrayUtil;
import org.metagene.genestrip.util.StringLongDigitTrie;
import org.metagene.genestrip.util.StringLongDigitTrie.StringLong;

public class KrakenResultFastqMerger {
	protected static final Log logger = LogFactory.getLog(KrakenResultFastqMerger.class);

	private final byte[] krakenChars;
	private final byte[] readDescriptor, outReadDescriptor;
	private final byte[] read;
	private final byte[] readProbs;

	private final BufferedLineReader bufferedLineReaderKraken;
	private final BufferedLineReader bufferedLineReaderFastQ;

	public KrakenResultFastqMerger(int maxReadSizeBytes) {
		krakenChars = new byte[maxReadSizeBytes];
		readDescriptor = new byte[maxReadSizeBytes];
		outReadDescriptor = new byte[maxReadSizeBytes];
		read = new byte[maxReadSizeBytes];
		readProbs = new byte[maxReadSizeBytes];

		bufferedLineReaderKraken = new BufferedLineReader();
		bufferedLineReaderFastQ = new BufferedLineReader();
	}

	public List<StringLong> process(InputStream fromKraken, InputStream fromFastQ,
			KrakenResultFastqMergeListener listener) throws IOException {
		StringLongDigitTrie root = new StringLongDigitTrie();

		int krakenPos;
		long readCount = 0;
		int pos;
		boolean start;
		boolean descriptor;
		boolean classId;
		boolean readSize;
		boolean fr;
		String classTaxid = null;
		int bps = 0;
		int startPos;
		int frStartPos;
		int readPos;

		bufferedLineReaderKraken.setInputStream(fromKraken);
		bufferedLineReaderFastQ.setInputStream(fromFastQ);

		for (krakenPos = bufferedLineReaderKraken
				.nextLine(krakenChars) - 1; krakenPos > 0; krakenPos = bufferedLineReaderKraken.nextLine(krakenChars) - 1) {
			krakenChars[krakenPos] = 0;

			if (fromFastQ != null) {
				pos = bufferedLineReaderFastQ.nextLine(readDescriptor) - 1;
				readDescriptor[pos] = 0;				
				
				pos = bufferedLineReaderFastQ.nextLine(read) - 1;
				read[pos] = 0;
				
				bufferedLineReaderFastQ.skipLine(); // Ignoring line 3.
				
				pos = bufferedLineReaderFastQ.nextLine(readProbs) - 1;
				readProbs[pos] = 0;
			}
			readCount++;

			start = true;
			descriptor = false;
			classId = false;
			readSize = false;
			fr = false;
			startPos = 0;
			frStartPos = 0;
			readPos = 0;

			for (int i = 0; i < krakenPos; i++) {
				if (krakenChars[i] == '\t') {
					if (start) {
						start = false;
						descriptor = true;
					} else if (descriptor) {
						descriptor = false;
						outReadDescriptor[i - 2] = 0;

						if (fromFastQ != null) {
							int j;
							for (j = 0; j < outReadDescriptor.length && outReadDescriptor[j] != 0; j++) {
								if (outReadDescriptor[j] != readDescriptor[j + 1]) {
//									System.out.println(ByteArrayUtil.toString(outReadDescriptor));
//									System.out.println(ByteArrayUtil.toString(readDescriptor));
									throw new IllegalStateException("Inconsistent files for read " + readCount);
								}
							}
							if (j == outReadDescriptor.length) {
								throw new IllegalStateException("Inconsistent kraken output...");
							}
						}

						classId = true;
						startPos = i + 1;
					} else if (classId) {
						classId = false;
						classTaxid = root.add(krakenChars, startPos, i, 0);
						readSize = true;
						startPos = i + 1;
					} else if (readSize) {
						readSize = false;
						bps = ByteArrayUtil.byteArrayToInt(krakenChars, startPos, i);
						startPos = i + 1;
					}
				} else if (descriptor) {
					outReadDescriptor[i - 2] = krakenChars[i];
				} else if (krakenChars[i] == ':') {
					fr = true;
					frStartPos = i + 1;
				} else if (fr && krakenChars[i] == ' ') {
					int frN = ByteArrayUtil.byteArrayToInt(krakenChars, frStartPos, i);
					String taxidStr = root.add(krakenChars, startPos, frStartPos - 1, frN);

					if (taxidStr != null && listener != null) {
						listener.newTaxIdForRead(readCount, readDescriptor, read, readProbs, classTaxid, bps, readPos,
								taxidStr, frN, krakenChars, root);
					}
					readPos += frN;

					startPos = i + 1;
				}
			}
			if (startPos < krakenPos && fr) {
				int frN = ByteArrayUtil.byteArrayToInt(krakenChars, frStartPos, krakenPos);
				String taxidStr = root.add(krakenChars, startPos, frStartPos - 1, frN);

				if (taxidStr != null && listener != null) {
					listener.newTaxIdForRead(readCount, readDescriptor, read, readProbs, classTaxid, bps, readPos,
							taxidStr, frN, krakenChars, root);
				}
			}
		}

		List<StringLong> list = new ArrayList<StringLongDigitTrie.StringLong>();
		root.collect(list);
		return list;
	}
}
