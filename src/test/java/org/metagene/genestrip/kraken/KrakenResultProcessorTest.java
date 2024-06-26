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

import static org.junit.Assert.assertEquals;

import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.io.InputStream;

import org.junit.Test;

public class KrakenResultProcessorTest {
	@Test
	public void testProcess() throws IOException {
		KrakenResultProcessor parser = new KrakenResultProcessor(4096);

		String testOut = "U	FP200005993L1C001R00807111253	0	150	0:89 3:27 11:2 0:18 1301:2 0:16 1301:3 0:5 1301:1 0:6 1301:5 0:53 28037:2 29606:1 0:5 9606:4 0:1 9606:3 0:20\n                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            "
				+ "C	A01245:102:H3JCKDMXY:1:1102:25789:122	9606	151	0:29 9606:5 0:6 9606:5 0:1 9606:2 0:8 9606:3 0:5 9606:1 0:12 9606:6 0:9 9606:1 0:24\n"
				+ "U	A01245:102:H3JCKDMXY:1:1102:23782:19413	0	151	0:23 28188:5 0:89\n"
				+ "C	A01245:102:H3JCKDMXY:1:1102:28664:19413	9606	151	0:5 9606:3 0:15 9606:1 0:2 9606:5 0:19 9606:1 0:12 9606:5 0:10 9606:1 0:31 9606:1 0:6\n"
				+ "C	FP200005993L1C001R00806844745	28037	150	0:1 1301:2 0:18 1301:2 0:16 1301:3 0:5 1301:1 0:6 1301:5 0:53 28037:2 29606:1 0:5 9606:4 0:1 9606:3 0:20\n";

		StringBuffer buffer = new StringBuffer();
		try (InputStream stream = new ByteArrayInputStream(testOut.getBytes(), 0, testOut.length())) {
			parser.process(stream, new KrakenResultListener() {
				@Override
				public void newTaxIdForRead(long lineCount, byte[] readDescriptor, String krakenTaxid, int bps, int pos,
						String kmerTaxid, int hitLength, byte[] output) {
					buffer.append(krakenTaxid);
					buffer.append(' ');
					buffer.append(bps);
					buffer.append(' ');
					buffer.append(pos);
					buffer.append(' ');
					buffer.append(kmerTaxid);
					buffer.append(' ');
					buffer.append(hitLength);
					buffer.append('\n');
				}
			});
		}

		String res = "0 150 0 0 89\n" + "0 150 89 3 27\n" + "0 150 116 11 2\n" + "0 150 118 0 18\n"
				+ "0 150 136 1301 2\n" + "0 150 138 0 16\n" + "0 150 154 1301 3\n" + "0 150 157 0 5\n"
				+ "0 150 162 1301 1\n" + "0 150 163 0 6\n" + "0 150 169 1301 5\n" + "0 150 174 0 53\n"
				+ "0 150 227 28037 2\n" + "0 150 229 29606 1\n" + "0 150 230 0 5\n" + "0 150 235 9606 4\n"
				+ "0 150 239 0 1\n" + "0 150 240 9606 3\n" + "0 150 243 0 20\n" + "9606 151 0 0 29\n"
				+ "9606 151 29 9606 5\n" + "9606 151 34 0 6\n" + "9606 151 40 9606 5\n" + "9606 151 45 0 1\n"
				+ "9606 151 46 9606 2\n" + "9606 151 48 0 8\n" + "9606 151 56 9606 3\n" + "9606 151 59 0 5\n"
				+ "9606 151 64 9606 1\n" + "9606 151 65 0 12\n" + "9606 151 77 9606 6\n" + "9606 151 83 0 9\n"
				+ "9606 151 92 9606 1\n" + "9606 151 93 0 24\n" + "0 151 0 0 23\n" + "0 151 23 28188 5\n"
				+ "0 151 28 0 89\n" + "9606 151 0 0 5\n" + "9606 151 5 9606 3\n" + "9606 151 8 0 15\n"
				+ "9606 151 23 9606 1\n" + "9606 151 24 0 2\n" + "9606 151 26 9606 5\n" + "9606 151 31 0 19\n"
				+ "9606 151 50 9606 1\n" + "9606 151 51 0 12\n" + "9606 151 63 9606 5\n" + "9606 151 68 0 10\n"
				+ "9606 151 78 9606 1\n" + "9606 151 79 0 31\n" + "9606 151 110 9606 1\n" + "9606 151 111 0 6\n"
				+ "28037 150 0 0 1\n" + "28037 150 1 1301 2\n" + "28037 150 3 0 18\n" + "28037 150 21 1301 2\n"
				+ "28037 150 23 0 16\n" + "28037 150 39 1301 3\n" + "28037 150 42 0 5\n" + "28037 150 47 1301 1\n"
				+ "28037 150 48 0 6\n" + "28037 150 54 1301 5\n" + "28037 150 59 0 53\n" + "28037 150 112 28037 2\n"
				+ "28037 150 114 29606 1\n" + "28037 150 115 0 5\n" + "28037 150 120 9606 4\n" + "28037 150 124 0 1\n"
				+ "28037 150 125 9606 3\n" + "28037 150 128 0 20\n";

//		System.out.println(buffer);

		assertEquals(res, buffer.toString());
	}
}
