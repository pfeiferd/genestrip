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
				+ "C	A01245:102:H3JCKDMXY:1:1102:28664:19413	9606	151	0:5 9606:3 0:15 9606:1 0:2 9606:5 0:19 9606:1 0:12 9606:5 0:10 9606:1 0:31 9606:1 0:6";

		InputStream stream = new ByteArrayInputStream(testOut.getBytes(StandardCharsets.UTF_8));

		parser.process(stream, new KrakenResultListener() {

			@Override
			public void newTaxIdForRead(long lineCount, byte[] readDescriptor, String krakenTaxid, int bps,
					String kmerTaxid, int hitLength) {
				System.out.print(krakenTaxid);
				System.out.print(' ');
				System.out.print(bps);
				System.out.print(' ');
				System.out.print(kmerTaxid);
				System.out.print(' ');
				System.out.println(hitLength);
			}
		});
	}
}
