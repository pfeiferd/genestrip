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

import java.io.PrintStream;
import java.util.HashSet;
import java.util.Set;

import org.metagene.genestrip.tax.TaxTree.TaxIdNode;
import org.metagene.genestrip.trie.KMerTrie;
import org.metagene.genestrip.util.ByteArrayToString;

public interface MergeListener {
	MergeListener STD_ERR_DEFAULT_UPDATE = new MergeListener() {
		@Override
		public void newTaxidForRead(long readCount, String taxid, byte[] readDescriptor, byte[] read,
				byte[] readProbs) {
			System.err.print("Read: ");
			System.err.println(readCount);
			System.err.print("New taxid returned: ");
			System.err.println(taxid);
			System.err.println(readDescriptor);
			System.err.println(read);
			System.err.println(readProbs);
		}
	};

	public void newTaxidForRead(long readCount, String taxid, byte[] readDescriptor, byte[] read, byte[] readProbs);

	public static MergeListener fillKMerTrie(KMerTrie<String> trie, MergeListener delegate) {
		return new MergeListener() {
			@Override
			public void newTaxidForRead(long readCount, String taxid, byte[] readDescriptor, byte[] read,
					byte[] readProbs) {
				trie.put(read, 0, taxid);
				if (readCount % 10000 == 0) {
					if (KrakenKMerFastqMerger.logger.isInfoEnabled()) {
						KrakenKMerFastqMerger.logger.info("Trie entries:" + trie.getEntries());
						KrakenKMerFastqMerger.logger.info("Trie put ratio:" + ((double) trie.getEntries() / readCount));
					}
				}
				if (delegate != null) {
					delegate.newTaxidForRead(readCount, taxid, readDescriptor, read, readProbs);
				}
			}
		};
	}

	public static MergeListener createFilterByTaxIdNodes(Set<TaxIdNode> taxIdNodes, MergeListener delegate) {
		Set<String> taxIds = new HashSet<String>();
		for (TaxIdNode node : taxIdNodes) {
			taxIds.add(node.getTaxId());
		}
		return createFilterByTaxId(taxIds, delegate);
	}

	public static MergeListener createFastQOutputFilterByTaxId(PrintStream printStream, MergeListener delegate) {
		return new MergeListener() {
			@Override
			public void newTaxidForRead(long readCount, String taxid, byte[] readDescriptor, byte[] read,
					byte[] readProbs) {
				ByteArrayToString.print(readDescriptor, printStream);
				printStream.print(":taxid:");
				printStream.println(taxid);
				ByteArrayToString.print(read, printStream);
				printStream.println();
				printStream.println("+");
				ByteArrayToString.print(readProbs, printStream);
				printStream.println();
	
				if (delegate != null) {
					delegate.newTaxidForRead(readCount, taxid, readDescriptor, read, readProbs);
				}
			}
		};
	}

	public static MergeListener createFilterByTaxId(final Set<String> taxIds, final MergeListener delegate) {
		return new MergeListener() {
			@Override
			public void newTaxidForRead(long readCount, String taxid, byte[] readDescriptor, byte[] read,
					byte[] readProbs) {
				if (taxIds.contains(taxid)) {
					delegate.newTaxidForRead(readCount, taxid, readDescriptor, read, readProbs);
				}
			}
		};
	}
}