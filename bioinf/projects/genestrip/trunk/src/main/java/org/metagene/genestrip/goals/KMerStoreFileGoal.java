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
import org.metagene.genestrip.kraken.KrakenResultFastqMergeListener;
import org.metagene.genestrip.kraken.KrakenResultFastqMerger;
import org.metagene.genestrip.make.FileGoal;
import org.metagene.genestrip.make.FileListGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.store.KMerStore;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;
import org.metagene.genestrip.util.ArraysUtil;
import org.metagene.genestrip.util.StreamProvider;

public class KMerStoreFileGoal extends FileListGoal<GSProject> {
	private final ObjectGoal<Set<TaxIdNode>, GSProject> taxNodesGoal;
	private final FileGoal<GSProject> krakenOutGoal;
	private final KMerFastqGoal kmerFastqGoal;
	private final ObjectGoal<Long, GSProject> sizeGoal;

	@SafeVarargs
	public KMerStoreFileGoal(GSProject project, String name, ObjectGoal<Set<TaxIdNode>, GSProject> taxNodesGoal,
			FileGoal<GSProject> krakenOutGoal, KMerFastqGoal kmerFastqGoal, ObjectGoal<Long, GSProject> sizeGoal,
			Goal<GSProject>... deps) {
		super(project, name, project.getOutputFile(name, FileType.SER),
				ArraysUtil.append(deps, taxNodesGoal, kmerFastqGoal, krakenOutGoal, sizeGoal));
		this.taxNodesGoal = taxNodesGoal;
		this.krakenOutGoal = krakenOutGoal;
		this.kmerFastqGoal = kmerFastqGoal;
		this.sizeGoal = sizeGoal;
	}

	@Override
	protected void makeFile(File trieFile) {
		try {
			KMerStore<String> store = getProject().getKMerStoreFactory().createKMerStore(String.class);
			store.initSize(sizeGoal.get());

			Set<String> taxIds = new HashSet<String>();
			for (TaxIdNode node : taxNodesGoal.get()) {
				taxIds.add(node.getTaxId());
			}

			KrakenResultFastqMergeListener filter = new KrakenResultFastqMergeListener() {
				private long counter;
				
				@Override
				public void newTaxIdForRead(long lineCount, byte[] readDescriptor, byte[] read, byte[] readProbs,
						String krakenTaxid, int bps, int pos, String kmerTaxid, int hitLength, byte[] output) {
					counter++;
					if (taxIds.contains(kmerTaxid)) {
						if (!store.put(read, 0, kmerTaxid, false)) {
							if (getLogger().isInfoEnabled()) {
								getLogger().info("Duplicate entry for read regarding taxid " + kmerTaxid);
							}
						}
						if (lineCount % 10000 == 0) {
							if (getLogger().isInfoEnabled()) {
								getLogger().info("Store entries:" + store.getEntries());
								getLogger().info("Store put ratio:" + ((double) store.getEntries() / counter));
							}
						}
					}
				}
			};

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
				getLogger().info("Stored entries: " + store.getEntries());
				getLogger().info("Saving File " + trieFile);
			}
			store.optimize();
			KMerStore.save(store, trieFile);
			if (getLogger().isInfoEnabled()) {
				getLogger().info("File saved " + trieFile);
			}
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}
}
