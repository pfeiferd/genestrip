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
package org.metagene.genestrip.fastq;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;

import java.io.IOException;

import org.junit.Test;
import org.metagene.genestrip.DefaultExecutionContext;
import org.metagene.genestrip.util.ByteArrayUtil;

public class FastqReaderTest {
	@Test
	public void testSingleThreadFastqReader() throws IOException {
		testSingleThreadFastqReaderProbs(true);
		testSingleThreadFastqReaderProbs(false);
	}

	protected void testSingleThreadFastqReaderProbs(boolean withProbs) throws IOException {
		int[] calls = new int[1];
		AbstractFastqReader fastqReader = new AbstractFastqReader(2, 3, 0, new DefaultExecutionContext(0, 1),
				withProbs) {
			@Override
			protected void nextEntry(ReadEntry readStruct, int threadIndex) throws IOException {
				calls[0]++;
				if (calls[0] == 1) {
					assertEquals("@S", ByteArrayUtil.toString(readStruct.readDescriptor));
					String s = "GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT";
					assertEquals(s, ByteArrayUtil.toString(readStruct.read));
					assertEquals(s.length(), readStruct.readSize);
					if (withProbs) {
						assertEquals("!''*((((***+))%%%++)(%%%%).1***-+*''))**55CCF>>>>>>CCCCCCC65",
								ByteArrayUtil.toString(readStruct.readProbs));
					} else {
						assertNull(readStruct.readProbs);
					}
					assertEquals(s.length(), readStruct.readProbsSize);
				} else if (calls[0] == 2) {
					assertEquals("@T", ByteArrayUtil.toString(readStruct.readDescriptor));
					String s = "CGAT";
					assertEquals(s, ByteArrayUtil.toString(readStruct.read));
					assertEquals(s.length(), readStruct.readSize);
					if (withProbs) {
						assertEquals("!**>", ByteArrayUtil.toString(readStruct.readProbs));
					} else {
						assertNull(readStruct.readProbs);
					}
					assertEquals(s.length(), readStruct.readProbsSize);
				}
			}
		};
		ClassLoader classLoader = getClass().getClassLoader();
		fastqReader.readFastq(classLoader.getResourceAsStream("fastq/SimpleTest.fastq"));
		assertEquals(2, calls[0]);
	}

	protected String probs(String s, boolean withProbs) {
		if (withProbs) {
			return s;
		}
		StringBuilder b = new StringBuilder();
		for (int i = 0; i < s.length(); i++) {
			b.append('~');
		}
		return b.toString();
	}
}
