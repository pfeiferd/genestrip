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
package org.metagene.genestrip.fastqgen;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;

import org.metagene.genestrip.io.StreamProvider;
import org.metagene.genestrip.store.KMerSortedArray;
import org.metagene.genestrip.store.KMerSortedArray.KMerSortedArrayVisitor;
import org.metagene.genestrip.util.CGAT;

public class KMerFastqGenerator {
	public static final String GENESTRIP_ID = "@GENESTRIP";

	private final KMerSortedArray<String> kmerStore;
	private final byte[] data;

	public KMerFastqGenerator(KMerSortedArray<String> kmerStore) {
		this.kmerStore = kmerStore;
		this.data = new byte[kmerStore.getK()];
	}

	public void generateFastq(File file, String taxid, String header) throws IOException {
		try (OutputStream out = StreamProvider.getOutputStreamForFile(file)) {
			generateFastq(out, taxid, header);
		}
	}
	
	public void generateFastq(OutputStream out, String taxid, String header) throws IOException {
		FastQWriter fastQWriter = new FastQWriter(GENESTRIP_ID + ":" + header, out);
		fastQWriter.start();
		kmerStore.visit(new KMerSortedArrayVisitor<String>() {
			@Override
			public void nextValue(KMerSortedArray<String> trie, long kmer, int index, long i) {
				String idFromDB = kmerStore.getValueForIndex(index);
				if (idFromDB.equals(taxid)) {
					fillData(kmer);
					fastQWriter.addRead(taxid, data);
				}
			}
		});
		fastQWriter.done();
	}
	
	protected void fillData(long kmer) {
		CGAT.longToKMerStraight(kmer, data, 0, data.length);
	}
}
