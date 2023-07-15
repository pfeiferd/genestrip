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
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.GSProject.FileType;
import org.metagene.genestrip.kraken.KrakenResultFastqMergeListener;
import org.metagene.genestrip.kraken.KrakenResultFastqMerger;
import org.metagene.genestrip.make.FileGoal;
import org.metagene.genestrip.make.FileListGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.store.KMerSortedArray;
import org.metagene.genestrip.store.KMerStoreWrapper;
import org.metagene.genestrip.store.KMerStoreWrapper.StoreStatsPerTaxid;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;
import org.metagene.genestrip.util.ArraysUtil;
import org.metagene.genestrip.util.DigitTrie;
import org.metagene.genestrip.util.StreamProvider;
import org.metagene.genestrip.util.StringLongDigitTrie;
import org.metagene.genestrip.util.StringLongDigitTrie.StringLong;

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
	protected void makeFile(File storeFile) {
		try {
			KMerSortedArray<String> store = (KMerSortedArray<String>) getProject().getKMerStoreFactory()
					.createKMerStore(String.class);
			store.initSize(sizeGoal.get());

			Set<String> taxIds = new HashSet<String>();
			for (TaxIdNode node : taxNodesGoal.get()) {
				taxIds.add(node.getTaxId());
			}

			System.out.println("***taxids***");
			System.out.println(taxIds);

			DigitTrie<StoreStatsPerTaxid> storeStats = new DigitTrie<StoreStatsPerTaxid>() {
				protected StoreStatsPerTaxid createInGet(String digits) {
					return new StoreStatsPerTaxid(digits);
				}
			};

			MyKrakenResultFastqMergeListener filter = new MyKrakenResultFastqMergeListener() {
				@Override
				public void newTaxIdForRead(long lineCount, byte[] readDescriptor, byte[] read, byte[] readProbs,
						String krakenTaxid, int bps, int pos, String kmerTaxid, int hitLength, byte[] output,
						StringLongDigitTrie root) {
					counter++;
					int colonPos = 0;
					int times = 0;
					for (colonPos = 0; readDescriptor[colonPos] != 0; colonPos++) {
						if (readDescriptor[colonPos] == ':') {
							times++;
							if (times == 3) {
								break;
							}
						}
					}
					if ("47466".equals(kmerTaxid)) {
						System.out.println("---found---");
					}
					String descriptorTaxid = null;
					if (times == 3) {
						int end;
						for (end = colonPos; readDescriptor[end] != 0; end++)
							;
						if ("47466".equals(kmerTaxid)) {
							System.out.println("---found again---");
							System.out.println(new String(readDescriptor, colonPos + 1, end -  colonPos + 1));
						}
						StringLong sl = root.get(readDescriptor, colonPos + 1, end, true);
						descriptorTaxid = sl.getStringValue();
						if (descriptorTaxid != null && taxIds.contains(descriptorTaxid)) {
							updateTotalKMers(descriptorTaxid);
						} else {
							if (getLogger().isInfoEnabled()) {
								getLogger().info("Bad taxid in read descriptor " + descriptorTaxid);
							}
						}
					} else {
						if (getLogger().isInfoEnabled()) {
							getLogger().info("Missing taxid in read descriptor");
						}
					}
					if (taxIds.contains(kmerTaxid)) {
						if (lastKMerTaxid == kmerTaxid) {
							contig++;
						} else {
							updateContigStats();
							contig = 0;
						}
						lastKMerTaxid = kmerTaxid;
						if (!store.put(read, 0, kmerTaxid, false)) {
							if (getLogger().isInfoEnabled()) {
								getLogger().info("Potential duplicate entry for kmer regarding taxid " + kmerTaxid);
							}
						} else {
							updateStoredKMers(descriptorTaxid);
						}
						if (lineCount % 1000000 == 0) {
							if (getLogger().isInfoEnabled()) {
								getLogger().info("Store entries: " + store.getEntries());
								getLogger().info("Store entry ratio: " + ((double) store.getEntries() / counter));
							}
						}
					} else {
						updateContigStats();
						contig = 0;
						lastKMerTaxid = null;
					}
					lastDescriptorTaxid = descriptorTaxid;
				}
			};
			filter.storeStats = storeStats;
			filter.totalStats = new StoreStatsPerTaxid(null);

			KrakenResultFastqMerger krakenKMerFastqMerger = new KrakenResultFastqMerger(
					getProject().getConfig().getMaxReadSizeBytes());

			for (int i = 0; i < krakenOutGoal.getFiles().size(); i++) {
				if (getLogger().isInfoEnabled()) {
					getLogger().info("Reading file " + krakenOutGoal.getFiles().get(i));
					getLogger().info("Reading file " + kmerFastqGoal.getFiles().get(i));
				}

				InputStream stream1 = StreamProvider.getInputStreamForFile(krakenOutGoal.getFiles().get(i));
				InputStream stream2 = StreamProvider.getInputStreamForFile(kmerFastqGoal.getFiles().get(i));

				filter.contig = 0;
				filter.lastKMerTaxid = null;
				filter.lastDescriptorTaxid = null;
				krakenKMerFastqMerger.process(stream1, stream2, filter);
				filter.updateContigStats();

				stream1.close();
				stream2.close();
			}

			if (getLogger().isInfoEnabled()) {
				getLogger().info("Total stored entries: " + store.getEntries());
				getLogger().info("Total store entry ratio: " + ((double) store.getEntries() / filter.counter));
				getLogger().info("Saving file " + storeFile);
			}
			store.optimize();

			List<StoreStatsPerTaxid> list = new ArrayList<StoreStatsPerTaxid>();
			storeStats.collect(list);
			list.add(filter.totalStats);

			KMerStoreWrapper wrapper = new KMerStoreWrapper((KMerSortedArray<String>) store, taxNodesGoal.get(), list);
			wrapper.save(storeFile);
			if (getLogger().isInfoEnabled()) {
				getLogger().info("File saved " + storeFile);
			}
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}

	private abstract static class MyKrakenResultFastqMergeListener implements KrakenResultFastqMergeListener {
		protected StoreStatsPerTaxid totalStats;
		protected DigitTrie<StoreStatsPerTaxid> storeStats;
		protected String lastKMerTaxid;
		protected long counter = 0;
		protected int contig;
		protected String lastDescriptorTaxid;

		public void updateContigStats() {
			if (contig > 0) {
				if (lastDescriptorTaxid != null) {
					StoreStatsPerTaxid stats = storeStats.get(lastDescriptorTaxid, true);
					stats.contigs++;
					if (contig > stats.maxContigLen) {
						stats.maxContigLen = contig;
					}
				}
			}
		}

		public void updateStoredKMers(String descriptorTaxid) {
			if (descriptorTaxid != null) {
				StoreStatsPerTaxid stats = storeStats.get(descriptorTaxid, true);
				stats.storedKMers++;
			}
			totalStats.storedKMers++;
		}

		public void updateTotalKMers(String descriptorTaxid) {
			StoreStatsPerTaxid stats = storeStats.get(descriptorTaxid, true);
			stats.totalKMers++;
			totalStats.totalKMers++;
		}
	}
}
