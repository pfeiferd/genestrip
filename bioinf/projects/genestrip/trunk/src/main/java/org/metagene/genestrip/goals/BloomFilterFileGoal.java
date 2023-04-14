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

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.util.HashSet;
import java.util.Set;

import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.bloom.KMerBloomIndex;
import org.metagene.genestrip.kraken.KrakenResultFastqMergeListener;
import org.metagene.genestrip.kraken.KrakenResultFastqMerger;
import org.metagene.genestrip.make.FileListGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;
import org.metagene.genestrip.util.ArraysUtil;
import org.metagene.genestrip.util.StreamProvider;

public class BloomFilterFileGoal extends FileListGoal<GSProject> {
	private final ObjectGoal<Set<TaxIdNode>, GSProject> taxNodesGoal;
	private final BloomFilterSizeGoal sizeGoal;

	@SafeVarargs
	public BloomFilterFileGoal(GSProject project, String name, BloomFilterSizeGoal sizeGoal,
			ObjectGoal<Set<TaxIdNode>, GSProject> taxNodesGoal, Goal<GSProject>... deps) {
		super(project, name, project.getBloomFilterFile(), ArraysUtil.append(deps, sizeGoal, taxNodesGoal));
		this.taxNodesGoal = taxNodesGoal;
		this.sizeGoal = sizeGoal;
	}

	@Override
	protected void makeFile(File bloomFilterFile) {
		try {
			KMerBloomIndex bloomIndex = new KMerBloomIndex(bloomFilterFile.getName(), getProject().getkMserSize(),
					sizeGoal.get(), 0.00001, null);

			if (getLogger().isInfoEnabled()) {
				getLogger().info("Number of k-mers for " + bloomIndex + ": " + sizeGoal.get());
				getLogger().info(
						"Bloom filter array size of " + bloomIndex + ": " + bloomIndex.getBitSize() / 8 / 1024 + "KB");
			}

			Set<String> taxIds = new HashSet<String>();
			for (TaxIdNode node : taxNodesGoal.get()) {
				taxIds.add(node.getTaxId());
			}

			KrakenResultFastqMergeListener filter = KrakenResultFastqMergeListener.createFilterByTaxIds(taxIds,
					new KrakenResultFastqMergeListener() {
						private long counter = 0;

						@Override
						public void newTaxIdForRead(long lineCount, byte[] readDescriptor, byte[] read,
								byte[] readProbs, String krakenTaxid, int bps, int pos, String kmerTaxid, int hitLength,
								byte[] output) {
							bloomIndex.putDirectKMer(read, 0);
							if (++counter % 10000 == 0) {
								if (getLogger().isInfoEnabled()) {
									getLogger().info("Added kmers to bloom filter: " + counter);
								}
							}
						}
					});
			KrakenResultFastqMerger krakenKMerFastqMerger = new KrakenResultFastqMerger(
					getProject().getConfig().getMaxReadSizeBytes());

			if (getLogger().isInfoEnabled()) {
				getLogger().info("Reading file " + getProject().getKrakenOutFile());
				getLogger().info("Reading file " + getProject().getKmerFastqFile());
			}

			InputStream stream1 = StreamProvider.getInputStreamForFile(getProject().getKrakenOutFile());
			InputStream stream2 = StreamProvider.getInputStreamForFile(getProject().getKmerFastqFile());
			krakenKMerFastqMerger.process(stream1, stream2, filter);
			stream1.close();
			stream2.close();

			if (getLogger().isInfoEnabled()) {
				getLogger().info("File save " + bloomFilterFile);
			}
			bloomIndex.save(bloomFilterFile);
			if (getLogger().isInfoEnabled()) {
				getLogger().info("File saved " + bloomFilterFile);
			}
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}

}
