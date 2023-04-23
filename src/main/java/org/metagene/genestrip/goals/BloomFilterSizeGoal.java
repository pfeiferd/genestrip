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
package org.metagene.genestrip.goals;

import java.io.IOException;
import java.io.InputStream;
import java.util.HashSet;
import java.util.Set;

import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.kraken.KrakenResultFastqMergeListener;
import org.metagene.genestrip.kraken.KrakenResultFastqMerger;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;
import org.metagene.genestrip.util.ArraysUtil;
import org.metagene.genestrip.util.StreamProvider;

public class BloomFilterSizeGoal extends ObjectGoal<Long, GSProject> {
	private final ObjectGoal<Set<TaxIdNode>, GSProject> taxNodesGoal;
	private final KrakenOutGoal krakenOutGoal;

	@SafeVarargs
	public BloomFilterSizeGoal(GSProject project, String name, ObjectGoal<Set<TaxIdNode>, GSProject> taxNodesGoal,
			KrakenOutGoal krakenOutGoal, Goal<GSProject>... deps) {
		super(project, name, ArraysUtil.append(deps, taxNodesGoal, krakenOutGoal));
		this.taxNodesGoal = taxNodesGoal;
		this.krakenOutGoal = krakenOutGoal;
	}

	@Override
	public void makeThis() {
		try {
			long[] counter = new long[] { 0 };
			Set<String> taxIds = new HashSet<String>();
			for (TaxIdNode node : taxNodesGoal.get()) {
				taxIds.add(node.getTaxId());
			}

			KrakenResultFastqMergeListener filter = KrakenResultFastqMergeListener.createFilterByTaxIds(taxIds,
					new KrakenResultFastqMergeListener() {
						@Override
						public void newTaxIdForRead(long lineCount, byte[] readDescriptor, byte[] read,
								byte[] readProbs, String krakenTaxid, int bps, int pos, String kmerTaxid, int hitLength,
								byte[] output) {
							counter[0]++;
						}
					});
			KrakenResultFastqMerger krakenKMerFastqMerger = new KrakenResultFastqMerger(
					getProject().getConfig().getMaxReadSizeBytes());

			if (getLogger().isInfoEnabled()) {
				getLogger().info("Reading file " + krakenOutGoal.getOutputFile());
			}
			InputStream stream = StreamProvider.getInputStreamForFile(krakenOutGoal.getOutputFile());
			krakenKMerFastqMerger.process(stream, null, filter);
			stream.close();

			set(counter[0]);
			if (getLogger().isInfoEnabled()) {
				getLogger().info("Count for kmer kraken out with taxids: " + counter[0]);
			}
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}
}
