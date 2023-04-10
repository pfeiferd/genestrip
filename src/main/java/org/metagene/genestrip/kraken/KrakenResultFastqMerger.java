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
import java.util.HashMap;
import java.util.Map;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.metagene.genestrip.util.BufferedLineReader;
import org.metagene.genestrip.util.CountingDigitTrie;

public class KrakenResultFastqMerger {
	protected static final Log logger = LogFactory.getLog(KrakenResultFastqMerger.class);

	private final byte[] krakenChars;
	private final byte[] readDescriptor, outReadDescriptor;
	private final byte[] read;
	private final byte[] readProbs;

	private final BufferedLineReader bufferedLineReaderKraken;
	private final BufferedLineReader bufferedLineReaderFastQ;

	public KrakenResultFastqMerger(int maxReadSizeBytes) {
		krakenChars = new byte[4096];
		readDescriptor = new byte[2048];
		outReadDescriptor = new byte[2048];
		read = new byte[maxReadSizeBytes];
		readProbs = new byte[maxReadSizeBytes];

		bufferedLineReaderKraken = new BufferedLineReader();
		bufferedLineReaderFastQ = new BufferedLineReader();
	}

	public Map<String, Long> process(InputStream bufferedInFromKraken, InputStream bufferedInFastQ,
			KrakenResultFastqMergeListener listener) throws IOException {
		CountingDigitTrie root = new CountingDigitTrie();

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

		bufferedLineReaderKraken.setInputStream(bufferedInFromKraken);
		bufferedLineReaderFastQ.setInputStream(bufferedInFastQ);

		for (krakenPos = bufferedLineReaderKraken
				.nextLine(krakenChars); krakenPos > 0; krakenPos = bufferedLineReaderKraken.nextLine(krakenChars)) {
			krakenChars[krakenPos - 1] = 0;

			if (bufferedInFastQ != null) {
				pos = bufferedLineReaderFastQ.nextLine(readDescriptor);
				readDescriptor[pos - 1] = 0;
				pos = bufferedLineReaderFastQ.nextLine(read);
				read[pos - 1] = 0;
				bufferedLineReaderFastQ.skipLine(); // Ignoring line 3.
				pos = bufferedLineReaderFastQ.nextLine(readProbs);
				readProbs[pos - 1] = 0;
			}

			start = true;
			descriptor = false;
			classId = false;
			readSize = false;
			fr = false;
			startPos = 0;
			frStartPos = 0;
			readPos = 0;

			for (int i = 0; i < krakenPos - 1; i++) {
				if (krakenChars[i] == '\t') {
					if (start) {
						start = false;
						descriptor = true;
					} else if (descriptor) {
						descriptor = false;
						outReadDescriptor[i - 2] = 0;

						if (bufferedInFastQ != null) {
							int j;
							for (j = 0; j < outReadDescriptor.length && outReadDescriptor[j] != 0; j++) {
								if (outReadDescriptor[j] != readDescriptor[j + 1]) {
									throw new IllegalStateException("In consistent files for read " + readCount);
								}
							}
							if (j == outReadDescriptor.length) {
								throw new IllegalStateException("In consistent kraken output...");
							}
						}

						classId = true;
						startPos = i + 1;
					} else if (classId) {
						classId = false;
						classTaxid = root.get(krakenChars, startPos, i - 1);
						readSize = true;
						startPos = i + 1;
					} else if (readSize) {
						readSize = false;
						bps = byteArrayToInt(krakenChars, startPos, i);
						startPos = i + 1;
					}
				} else if (descriptor) {
					outReadDescriptor[i - 2] = krakenChars[i];
				} else if (krakenChars[i] == ':') {
					fr = true;
					frStartPos = i + 1;
				} else if (fr && krakenChars[i] == ' ') {
					int frN = byteArrayToInt(krakenChars, frStartPos, i);
					String taxidStr = root.add(krakenChars, startPos, frStartPos - 2, frN);

					if (listener != null) {
						listener.newTaxIdForRead(readCount, readDescriptor, read, readProbs, classTaxid, bps, readPos,
								taxidStr, frN, krakenChars);
					}
					readPos += frN;

					startPos = i + 1;
				}
			}
			if (startPos < krakenPos - 1 && fr) {
				int frN = byteArrayToInt(krakenChars, frStartPos, krakenPos - 1);
				String taxidStr = root.add(krakenChars, startPos, frStartPos - 2, frN);

				if (listener != null) {
					listener.newTaxIdForRead(readCount, readDescriptor, read, readProbs, classTaxid, bps, readPos,
							taxidStr, frN, krakenChars);
				}
			}
			readCount++;
		}

		Map<String, Long> map = new HashMap<String, Long>();
		root.collect(map);
		return map;
	}

	private int byteArrayToInt(byte[] data, int start, int end) throws NumberFormatException {
		int result = 0;
		for (int i = start; i < end; i++) {
			int digit = data[i] - '0';
			if ((digit < 0) || (digit > 9)) {
				System.out.println(new String(data));
				throw new NumberFormatException();
			}
			result *= 10;
			result += digit;
		}
		return result;
	}
}
