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
package org.metagene.genestrip.goals.refseq;

import java.io.File;
import java.io.IOException;
import java.util.Map;
import java.util.Set;

import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.GSProject.FileType;
import org.metagene.genestrip.bloom.MurmurCGATBloomFilter;
import org.metagene.genestrip.make.FileListGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.refseq.AbstractStoreFastaReader;
import org.metagene.genestrip.refseq.AccessionMap;
import org.metagene.genestrip.refseq.RefSeqCategory;
import org.metagene.genestrip.store.KMerSortedArray;
import org.metagene.genestrip.store.KMerStoreWrapper;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;

public class FillStoreGoal extends FileListGoal<GSProject> {
	private final ObjectGoal<Set<RefSeqCategory>[], GSProject> categoriesGoal;
	private final ObjectGoal<Set<TaxIdNode>, GSProject> taxNodesGoal;
	private final RefSeqFnaFilesDownloadGoal fnaFilesGoal;
	private final ObjectGoal<Map<File, TaxIdNode>, GSProject> additionalGoal;
	private final ObjectGoal<AccessionMap, GSProject> accessionMapGoal;
	private final ObjectGoal<MurmurCGATBloomFilter, GSProject> bloomFilterGoal;
	private FilledStoreGoal filledStoreGoal;

	@SafeVarargs
	public FillStoreGoal(GSProject project, String name, ObjectGoal<Set<RefSeqCategory>[], GSProject> categoriesGoal,
			ObjectGoal<Set<TaxIdNode>, GSProject> taxNodesGoal, RefSeqFnaFilesDownloadGoal fnaFilesGoal,
			ObjectGoal<Map<File, TaxIdNode>, GSProject> additionalGoal,
			ObjectGoal<AccessionMap, GSProject> accessionMapGoal, FillSizeGoal fillSizeGoal,
			ObjectGoal<MurmurCGATBloomFilter, GSProject> bloomFilterGoal, Goal<GSProject>... deps) {
		super(project, name, project.getOutputFile(name, FileType.DB), Goal.append(deps, categoriesGoal, taxNodesGoal,
				fnaFilesGoal, accessionMapGoal, fillSizeGoal, bloomFilterGoal));
		this.categoriesGoal = categoriesGoal;
		this.taxNodesGoal = taxNodesGoal;
		this.fnaFilesGoal = fnaFilesGoal;
		this.additionalGoal = additionalGoal;
		this.accessionMapGoal = accessionMapGoal;
		this.bloomFilterGoal = bloomFilterGoal;
	}

	public void setFilledStoreGoal(FilledStoreGoal filledStoreGoal) {
		this.filledStoreGoal = filledStoreGoal;
	}

	@Override
	public void makeFile(File storeFile) {
		KMerSortedArray<String> store = new KMerSortedArray<String>(getProject().getConfig().getKMerSize(),
				getProject().getConfig().getBloomFilterFpp(), null, false);
		store.initSize(bloomFilterGoal.get().getEntries());

		try {
			MyFastaReader fastaReader = new MyFastaReader(getProject().getConfig().getMaxReadSizeBytes(),
					taxNodesGoal.get(), accessionMapGoal.get(), store, getProject().getMaxGenomesPerTaxid(), getProject().getMaxDust());

			for (File fnaFile : fnaFilesGoal.getFiles()) {
				RefSeqCategory cat = fnaFilesGoal.getCategoryForFile(fnaFile);
				if (categoriesGoal.get()[0].contains(cat)) {
					fastaReader.readFasta(fnaFile);
				}
			}
			Map<File, TaxIdNode> additionalMap = additionalGoal.get();
			for (File additionalFasta : additionalMap.keySet()) {
				TaxIdNode node = additionalMap.get(additionalFasta);
				if (taxNodesGoal.get().contains(node)) {
					fastaReader.ignoreAccessionMap(node);
					fastaReader.readFasta(additionalFasta);
				}
			}
			if (getLogger().isWarnEnabled()) {
				getLogger().warn("Not stored kmers: " + fastaReader.tooManyCounter);
			}

			store.optimize();

			KMerStoreWrapper wrapper = new KMerStoreWrapper((KMerSortedArray<String>) store);
			File filterFile = getProject().getFilterFile(this);
			wrapper.save(storeFile, filterFile);
			if (getLogger().isInfoEnabled()) {
				getLogger().info("File saved " + storeFile + " along with index " + filterFile);
			}
			filledStoreGoal.setStoreWrapper(wrapper);
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}

	protected static class MyFastaReader extends AbstractStoreFastaReader {
		private final KMerSortedArray<String> store;
		private long tooManyCounter;

		public MyFastaReader(int bufferSize, Set<TaxIdNode> taxNodes, AccessionMap accessionMap,
				KMerSortedArray<String> store, int maxGenomesPerTaxId, int maxDust) {
			super(bufferSize, taxNodes, accessionMap, store.getK(), maxGenomesPerTaxId, maxDust);
			this.store = store;
		}

		@Override
		protected void handleStore() {
			if (store.isFull()) {
				tooManyCounter++;
			} else {
				if (node.getTaxId() != null) {
					store.putLong(byteRingBuffer.getKMer() , node.getTaxId());
				}
				else {
					if (getLogger().isWarnEnabled()) {
						getLogger().warn("Tax id node without taxid: " + node.getName());
					}
				}
			}
		}
	}
}