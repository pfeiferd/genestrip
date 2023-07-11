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
import org.metagene.genestrip.GSProject.FileType;
import org.metagene.genestrip.bloom.AbstractKMerBloomIndex;
import org.metagene.genestrip.bloom.MurmurCGATBloomFilter;
import org.metagene.genestrip.kraken.KrakenResultFastqMergeListener;
import org.metagene.genestrip.kraken.KrakenResultFastqMerger;
import org.metagene.genestrip.make.FileGoal;
import org.metagene.genestrip.make.FileListGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;
import org.metagene.genestrip.util.ArraysUtil;
import org.metagene.genestrip.util.CountingDigitTrie;
import org.metagene.genestrip.util.StreamProvider;

public class BloomFilterFileGoal extends FileListGoal<GSProject> {
	private final ObjectGoal<Set<TaxIdNode>, GSProject> taxNodesGoal;
	private final BloomFilterSizeGoal sizeGoal;
	private final FileGoal<GSProject> krakenOutGoal;
	private final KMerFastqGoal kmerFastqGoal;

	@SafeVarargs
	public BloomFilterFileGoal(GSProject project, String name, BloomFilterSizeGoal sizeGoal,
			ObjectGoal<Set<TaxIdNode>, GSProject> taxNodesGoal, FileGoal<GSProject> krakenOutGoal,
			KMerFastqGoal kmerFastqGoal, Goal<GSProject>... deps) {
		super(project, name, project.getOutputFile(name, FileType.SER),
				ArraysUtil.append(deps, sizeGoal, taxNodesGoal, krakenOutGoal, kmerFastqGoal));
		this.taxNodesGoal = taxNodesGoal;
		this.sizeGoal = sizeGoal;
		this.krakenOutGoal = krakenOutGoal;
		this.kmerFastqGoal = kmerFastqGoal;
	}

	@Override
	protected void makeFile(File bloomFilterFile) {
		try {

			AbstractKMerBloomIndex bloomIndex = new AbstractKMerBloomIndex(bloomFilterFile.getName(),
					getProject().getKMserSize(), 0.0001, null);
			// I found this out by trial end error: Guava bloom filter cant keep the FP-rate
			// when total entry size is too small.
			// So keep it to a minimum.
			// Maybe there is more (bad stuff) to it.
			long size = Math.max(1000 * 10000, sizeGoal.get());
			bloomIndex.getFilter().clear();
			bloomIndex.getFilter().ensureExpectedSize(size, false);

			if (getLogger().isInfoEnabled()) {
				getLogger().info("Number of k-mers for " + bloomIndex + ": " + sizeGoal.get());
				getLogger().info("Bloom filter array size of " + bloomIndex + ": "
						+ bloomIndex.getFilter().getBitSize() / 8 /  1024 + "KB");
			}

			Set<String> taxIds = new HashSet<String>();
			for (TaxIdNode node : taxNodesGoal.get()) {
				taxIds.add(node.getTaxId());
			}

			MurmurCGATBloomFilter bloomFilter = bloomIndex.getFilter();

			KrakenResultFastqMergeListener filter = KrakenResultFastqMergeListener.createFilterByTaxIds(taxIds,
					new KrakenResultFastqMergeListener() {
						private long counter = 0;

						@Override
						public void newTaxIdForRead(long lineCount, byte[] readDescriptor, byte[] read,
								byte[] readProbs, String krakenTaxid, int bps, int pos, String kmerTaxid, int hitLength,
								byte[] output, CountingDigitTrie root) {
							bloomFilter.put(read, 0);
							if (++counter % 10000 == 0) {
								if (getLogger().isInfoEnabled()) {
									getLogger().info("Added kmers to bloom filter: " + counter);
								}
							}
						}
					});
			KrakenResultFastqMerger krakenKMerFastqMerger = new KrakenResultFastqMerger(
					getProject().getConfig().getMaxReadSizeBytes());

			for (int i = 0; i < krakenOutGoal.getFiles().size(); i++) {
				if (getLogger().isInfoEnabled()) {
					getLogger().info("Reading file " + krakenOutGoal.getFiles().get(i));
					getLogger().info("Reading file " + kmerFastqGoal.getFiles().get(i));
				}

				InputStream stream1 = StreamProvider.getInputStreamForFile(krakenOutGoal.getFiles().get(i));
				InputStream stream2 = StreamProvider.getInputStreamForFile(kmerFastqGoal.getFiles().get(i));
				krakenKMerFastqMerger.process(stream1, stream2, filter);
				stream1.close();
				stream2.close();
			}

			if (getLogger().isInfoEnabled()) {
				getLogger().info("Saving file " + bloomFilterFile);
			}
			bloomIndex.save(bloomFilterFile);
			if (getLogger().isInfoEnabled()) {
				getLogger().info("Saved file " + bloomFilterFile);
			}
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}

}
