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

import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.nio.charset.StandardCharsets;

import junit.framework.TestCase;

public class KrakenResultParserTest extends TestCase {
	public void testProcess() throws IOException {
		KrakenResultParser parser = new KrakenResultParser();

		String testOut = "C	A01245:102:H3JCKDMXY:1:1102:25789:122	9606	151	0:29 9606:5 0:6 9606:5 0:1 9606:2 0:8 9606:3 0:5 9606:1 0:12 9606:6 0:9 9606:1 0:24\n"
				+ "U	A01245:102:H3JCKDMXY:1:1102:23782:19413	0	151	0:23 28188:5 0:89\n"
				+ "C	A01245:102:H3JCKDMXY:1:1102:28664:19413	9606	151	0:5 9606:3 0:15 9606:1 0:2 9606:5 0:19 9606:1 0:12 9606:5 0:10 9606:1 0:31 9606:1 0:6\n"
				+ "C	FP200005993L1C001R00806844745	28037	150	0:1 1301:2 0:18 1301:2 0:16 1301:3 0:5 1301:1 0:6 1301:5 0:53 28037:2 0:29606:1 0:5 9606:4 0:1 9606:3 0:20";

		InputStream stream = new ByteArrayInputStream(testOut.getBytes(StandardCharsets.UTF_8));
		
		StringBuffer buffer = new StringBuffer();

		parser.process(stream, new KrakenResultListener() {

			@Override
			public void newTaxIdForRead(long lineCount, byte[] readDescriptor, String krakenTaxid, int bps,
					String kmerTaxid, int hitLength) {
				buffer.append(krakenTaxid);
				buffer.append(' ');
				buffer.append(bps);
				buffer.append(' ');
				buffer.append(kmerTaxid);
				buffer.append(' ');
				buffer.append(hitLength);
				buffer.append('\n');
			}
		});
		
		String res = "9606 151 0 29\n"
				+ "9606 151 9606 5\n"
				+ "9606 151 0 6\n"
				+ "9606 151 9606 5\n"
				+ "9606 151 0 1\n"
				+ "9606 151 9606 2\n"
				+ "9606 151 0 8\n"
				+ "9606 151 9606 3\n"
				+ "9606 151 0 5\n"
				+ "9606 151 9606 1\n"
				+ "9606 151 0 12\n"
				+ "9606 151 9606 6\n"
				+ "9606 151 0 9\n"
				+ "9606 151 9606 1\n"
				+ "9606 151 0 24\n"
				+ "0 151 0 23\n"
				+ "0 151 28188 5\n"
				+ "0 151 0 89\n"
				+ "9606 151 0 5\n"
				+ "9606 151 9606 3\n"
				+ "9606 151 0 15\n"
				+ "9606 151 9606 1\n"
				+ "9606 151 0 2\n"
				+ "9606 151 9606 5\n"
				+ "9606 151 0 19\n"
				+ "9606 151 9606 1\n"
				+ "9606 151 0 12\n"
				+ "9606 151 9606 5\n"
				+ "9606 151 0 10\n"
				+ "9606 151 9606 1\n"
				+ "9606 151 0 31\n"
				+ "9606 151 9606 1\n"
				+ "9606 151 0 6\n"
				+ "28037 150 0 1\n"
				+ "28037 150 1301 2\n"
				+ "28037 150 0 18\n"
				+ "28037 150 1301 2\n"
				+ "28037 150 0 16\n"
				+ "28037 150 1301 3\n"
				+ "28037 150 0 5\n"
				+ "28037 150 1301 1\n"
				+ "28037 150 0 6\n"
				+ "28037 150 1301 5\n"
				+ "28037 150 0 53\n"
				+ "28037 150 28037 2\n"
				+ "28037 150 0 5\n"
				+ "28037 150 9606 4\n"
				+ "28037 150 0 1\n"
				+ "28037 150 9606 3\n"
				+ "28037 150 0 20\n";
		
//		System.out.println(buffer);
		
		assertEquals(res, buffer.toString());
	}
}
