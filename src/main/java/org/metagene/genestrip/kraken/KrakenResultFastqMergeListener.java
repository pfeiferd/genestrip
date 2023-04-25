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
import org.metagene.genestrip.util.ByteArrayUtil;

public interface KrakenResultFastqMergeListener {
	public void newTaxIdForRead(long lineCount, byte[] readDescriptor, byte[] read, byte[] readProbs,
			String krakenTaxid, int bps, int pos, String kmerTaxid, int hitLength, byte[] output);

	public static KrakenResultFastqMergeListener createPrintListener(PrintStream out,
			KrakenResultFastqMergeListener delegate, boolean withKrakenOut) {
		return new KrakenResultFastqMergeListener() {
			@Override
			public void newTaxIdForRead(long lineCount, byte[] readDescriptor, byte[] read, byte[] readProbs,
					String krakenTaxid, int bps, int pos, String kmerTaxid, int hitLength, byte[] output) {
				ByteArrayUtil.print(readDescriptor, out);
				out.println();
				ByteArrayUtil.print(read, out);
				out.println();
				out.println("+");
				ByteArrayUtil.print(readProbs, out);
				out.println();
				if (withKrakenOut) {
					ByteArrayUtil.print(output, out);
					out.println();
				}
				if (delegate != null) {
					delegate.newTaxIdForRead(lineCount, readDescriptor, read, readProbs, krakenTaxid, bps, pos,
							kmerTaxid, hitLength, output);
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
		return createFilterByTaxIds(taxIds, delegate);
	}

	public static KrakenResultFastqMergeListener createFilterByTaxIds(final Set<String> taxIds,
			final KrakenResultFastqMergeListener delegate) {
		return new KrakenResultFastqMergeListener() {
			@Override
			public void newTaxIdForRead(long lineCount, byte[] readDescriptor, byte[] read, byte[] readProbs,
					String krakenTaxid, int bps, int pos, String kmerTaxid, int hitLength, byte[] output) {
				if (taxIds.contains(kmerTaxid)) {
					delegate.newTaxIdForRead(lineCount, readDescriptor, read, readProbs, krakenTaxid, bps, pos,
							kmerTaxid, hitLength, output);
				}
			}
		};
	}
}
