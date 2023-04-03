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

public interface KrakenResultFastqMergeListener {
	public void newTaxIdForRead(long lineCount, byte[] readDescriptor, byte[] read, byte[] readProbs,
			String krakenTaxid, int bps, String kmerTaxid, int hitLength, byte[] output);

	public static KrakenResultFastqMergeListener fillKMerTrie(KMerTrie<String> trie,
			KrakenResultFastqMergeListener delegate) {
		return new KrakenResultFastqMergeListener() {
			@Override
			public void newTaxIdForRead(long lineCount, byte[] readDescriptor, byte[] read, byte[] readProbs,
					String krakenTaxid, int bps, String kmerTaxid, int hitLength, byte[] output) {
				trie.put(read, 0, kmerTaxid);
				if (lineCount % 10000 == 0) {
					if (KrakenKMerFastqMerger.logger.isInfoEnabled()) {
						KrakenKMerFastqMerger.logger.info("Trie entries:" + trie.getEntries());
						KrakenKMerFastqMerger.logger.info("Trie put ratio:" + ((double) trie.getEntries() / lineCount));
					}
				}
				if (delegate != null) {
					delegate.newTaxIdForRead(lineCount, readDescriptor, read, readProbs, krakenTaxid, bps, kmerTaxid,
							hitLength, output);
				}
			}
		};
	}

	public static KrakenResultFastqMergeListener createFilterByTaxIdNodes(Set<TaxIdNode> taxIdNodes,
			KrakenResultFastqMergeListener delegate) {
		Set<String> taxIds = new HashSet<String>();
		for (TaxIdNode node : taxIdNodes) {
			taxIds.add(node.getTaxId());
		}
		return createFilterByTaxId(taxIds, delegate);
	}

	public static KrakenResultFastqMergeListener createFastQOutputFilterByTaxId(PrintStream printStream,
			KrakenResultFastqMergeListener delegate) {
		return new KrakenResultFastqMergeListener() {
			private long lastLineCount;
			
			@Override
			public void newTaxIdForRead(long lineCount, byte[] readDescriptor, byte[] read, byte[] readProbs,
					String krakenTaxid, int bps, String kmerTaxid, int hitLength, byte[] output) {
				if (lastLineCount == lineCount) {
					throw new IllegalStateException("too many k-mers for read");
				}
				lastLineCount = lineCount;
				ByteArrayToString.print(readDescriptor, printStream);
				printStream.print(":taxid:");
				printStream.println(kmerTaxid);
				ByteArrayToString.print(read, printStream);
				printStream.println();
				printStream.println("+");
				ByteArrayToString.print(readProbs, printStream);
				printStream.println();

				if (delegate != null) {
					delegate.newTaxIdForRead(lineCount, readDescriptor, read, readProbs, krakenTaxid, bps, kmerTaxid,
							hitLength, output);
					;
				}
			}
		};
	}

	public static KrakenResultFastqMergeListener createFilterByTaxId(final Set<String> taxIds,
			final KrakenResultFastqMergeListener delegate) {
		return new KrakenResultFastqMergeListener() {
			@Override
			public void newTaxIdForRead(long lineCount, byte[] readDescriptor, byte[] read, byte[] readProbs,
					String krakenTaxid, int bps, String kmerTaxid, int hitLength, byte[] output) {
				if (taxIds.contains(kmerTaxid)) {
					delegate.newTaxIdForRead(lineCount, readDescriptor, read, readProbs, krakenTaxid, bps, kmerTaxid,
							hitLength, output);
				}
			}
		};
	}
}
