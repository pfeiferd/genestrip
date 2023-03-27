package org.metagene.genestrip.kraken;

import java.io.IOException;
import java.io.InputStream;
import java.util.HashMap;
import java.util.Map;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.metagene.genestrip.util.CountingDigitTrie;

public class KrakenResultParser {
	protected static final Log logger = LogFactory.getLog(KrakenKMerFastqMerger.class);

	private final byte[] krakenChars;
	private final byte[] readDescriptor;

	public KrakenResultParser() {
		krakenChars = new byte[4096];
		readDescriptor = new byte[4096];
	}

	public Map<String, Long> process(InputStream bufferedInFromKraken, KrakenResultListener listener)
			throws IOException {
		CountingDigitTrie root = new CountingDigitTrie();

		int krakenPos;
		long readCount = 0;
		for (int c = bufferedInFromKraken.read(); c != -1; c = bufferedInFromKraken.read()) {
			for (krakenPos = 0; c != -1 && c != '\n'; krakenPos++) {
				krakenChars[krakenPos] = (byte) c;
				c = bufferedInFromKraken.read();
			}

			boolean start = true;
			boolean descriptor = false;
			boolean classId = false;
			boolean readSize = false;
			boolean fr = false;
			String classTaxid = null;
			int bps = 0;
			int startPos = 0;
			int frStartPos = 0;
			for (int i = 0; i < krakenPos; i++) {
				if (krakenChars[i] == '\t') {
					if (start) {
						start = false;
						descriptor = true;
					} else if (descriptor) {
						descriptor = false;
						readDescriptor[i] = 0;
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
					readDescriptor[i] = krakenChars[i];
				} else if (krakenChars[i] == ':') {
					fr = true;
					frStartPos = i + 1;
				} else if (fr && krakenChars[i] == ' ') {
					int frN = byteArrayToInt(krakenChars, frStartPos, i);
					String taxidStr = root.add(krakenChars, startPos, frStartPos - 2, frN);

					if (listener != null) {
						listener.newTaxIdForRead(readCount, readDescriptor, classTaxid, bps, taxidStr, frN);
					}

					startPos = i + 1;
				}
			}
			int frN = byteArrayToInt(krakenChars, frStartPos, krakenPos);
			String taxidStr = root.add(krakenChars, startPos, frStartPos - 2, frN);

			if (listener != null) {
				listener.newTaxIdForRead(readCount, readDescriptor, classTaxid, bps, taxidStr, frN);
			}

			readCount++;
		}
		bufferedInFromKraken.close();

		Map<String, Long> map = new HashMap<String, Long>();
		root.collect(map);
		return map;
	}

	private int byteArrayToInt(byte[] data, int start, int end) throws NumberFormatException {
		int result = 0;
		for (int i = start; i < end; i++) {
			int digit = data[i] - '0';
			if ((digit < 0) || (digit > 9))
				throw new NumberFormatException();
			result *= 10;
			result += digit;
		}
		return result;
	}
}
