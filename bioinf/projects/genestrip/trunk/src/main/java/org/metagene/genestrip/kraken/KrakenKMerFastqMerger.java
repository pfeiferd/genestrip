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
import org.metagene.genestrip.util.CountingDigitTrie;

public class KrakenKMerFastqMerger {
	protected static final Log logger = LogFactory.getLog(KrakenKMerFastqMerger.class);

	private final byte[] krakenChars;
	private final byte[] readDescriptor;
	private final byte[] read;
	private final byte[] readProbs;

	public KrakenKMerFastqMerger(int maxReadSizeBytes) {
		krakenChars = new byte[2048];
		readDescriptor = new byte[2048];
		read = new byte[maxReadSizeBytes];
		readProbs = new byte[maxReadSizeBytes];
	}

	public Map<String, Long> process(InputStream bufferedInFromKraken, InputStream bufferedInFastQ,
			MergeListener update) throws IOException {
		CountingDigitTrie root = new CountingDigitTrie();

		int krakenPos;
		long readCount = 0;
		for (int c = bufferedInFromKraken.read(), d = bufferedInFastQ == null ? 0 : bufferedInFastQ.read(); c != -1
				&& d != -1; c = bufferedInFromKraken.read(), d = bufferedInFastQ == null ? 0 : bufferedInFastQ.read()) {
			for (krakenPos = 0; c != -1 && c != '\n'; krakenPos++) {
				krakenChars[krakenPos] = (byte) c;
				c = bufferedInFromKraken.read();
			}

			if (bufferedInFastQ != null) {
				int pos;
				for (pos = 0; d != -1 && d != '\n'; pos++) {
					readDescriptor[pos] = (byte) d;
					d = bufferedInFastQ.read();
				}
				readDescriptor[pos] = 0;
				d = bufferedInFastQ.read();
				for (pos = 0; d != -1 && d != '\n'; pos++) {
					read[pos] = (byte) d;
					d = bufferedInFastQ.read();
				}
				read[pos] = 0;
				// Ignore line with "+...":
				d = bufferedInFastQ.read();
				for (; d != -1 && d != '\n';) {
					d = bufferedInFastQ.read();
				}
				d = bufferedInFastQ.read();
				for (pos = 0; d != -1 && d != '\n'; pos++) {
					readProbs[pos] = (byte) d;
					d = bufferedInFastQ.read();
				}
				readProbs[pos] = 0;
			}
			readCount++;
			krakenPos--;
			int start;
			int end = krakenPos - 2;
			for (start = krakenPos; start > 0; start--) {
				if (krakenChars[start] == ':') {
					end = krakenPos;
				}
				if (krakenChars[start] == '\t') {
					break;
				}
			}
			String taxid = root.inc(krakenChars, start + 1, end - 2);

			if (bufferedInFastQ != null) {
				int i;
				for (i = 1; i < krakenChars.length && krakenChars[i + 1] != '\t'; i++) {
					if (krakenChars[i + 1] != readDescriptor[i]) {
						throw new IllegalStateException("In consistent files for read " + readCount);
					}
				}
				if (i == krakenChars.length) {
					throw new IllegalStateException("In consistent kraken output...");
				}
			}

			if (taxid != null && update != null) {
				update.newTaxidForRead(readCount, taxid, readDescriptor, read, readProbs);
			}
		}
		if (bufferedInFastQ != null) {
			bufferedInFastQ.close();
		}
		bufferedInFromKraken.close();

		Map<String, Long> map = new HashMap<String, Long>();
		root.collect(map);
		return map;
	}
}
