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
import org.metagene.genestrip.store.Database;
import org.metagene.genestrip.store.KMerStore;
import org.metagene.genestrip.store.KMerStore.IndexedKMerStoreVisitor;
import org.metagene.genestrip.tax.SmallTaxTree;
import org.metagene.genestrip.tax.TaxTree;
import org.metagene.genestrip.util.CGAT;

public class KMerFastqGenerator {
	public static final String GENESTRIP_ID = "@GENESTRIP";

	// private final Database database;
	private final KMerStore<SmallTaxTree.SmallTaxIdNode> kmerStore;
	private final byte[] data;

	public KMerFastqGenerator(Database database) {
		//this.database = database;
		this.data = new byte[database.getKmerStore().getK()];
		kmerStore = database.convertKMerStore();
	}

	public void generateFastq(File file, String taxid, String header, boolean withDesc) throws IOException {
		try (OutputStream out = StreamProvider.getOutputStreamForFile(file)) {
			generateFastq(out, taxid, header, withDesc);
		}
	}
	
	public void generateFastq(OutputStream out, String taxid, String header, boolean withDesc) throws IOException {
		FastQWriter fastQWriter = new FastQWriter(GENESTRIP_ID + ":" + header, out);
		fastQWriter.start();

		kmerStore.visit(new IndexedKMerStoreVisitor<SmallTaxTree.SmallTaxIdNode>() {
			@Override
			public void nextValue(KMerStore<SmallTaxTree.SmallTaxIdNode> trie, long kmer, int index, long i) {
				SmallTaxTree.SmallTaxIdNode node = kmerStore.getValueForIndex(index);
				if (isMatchingNode(node, taxid, withDesc)) {
					CGAT.longToKMerStraight(kmer, data, 0, data.length);
					fastQWriter.addRead(node.getTaxId(), data);
				}
			}
		});
		fastQWriter.done();
	}

	protected boolean isMatchingNode(SmallTaxTree.SmallTaxIdNode node, String taxid, boolean withDesc) {
		if (taxid == null) {
			return withDesc;
		}
		else if (withDesc) {
			while (node != null) {
				if (node.getTaxId().equals(taxid)) {
					return true;
				}
				node = node.getParent();
			}
			return false;
		}
		else {
			return node.getTaxId().equals(taxid);
		}
	}
}
