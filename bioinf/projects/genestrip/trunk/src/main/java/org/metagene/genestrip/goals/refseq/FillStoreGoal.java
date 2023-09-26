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
import java.util.Set;

import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.GSProject.FileType;
import org.metagene.genestrip.bloom.MurmurCGATBloomFilter;
import org.metagene.genestrip.goals.AdditionalFastasGoal;
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
	private final AdditionalFastasGoal additionalGoal;
	private final ObjectGoal<AccessionMap, GSProject> accessionMapGoal;
	private final ObjectGoal<MurmurCGATBloomFilter, GSProject> bloomFilterGoal;
	private FilledStoreGoal filledStoreGoal;

	@SafeVarargs
	public FillStoreGoal(GSProject project, String name, ObjectGoal<Set<RefSeqCategory>[], GSProject> categoriesGoal,
			ObjectGoal<Set<TaxIdNode>, GSProject> taxNodesGoal, RefSeqFnaFilesDownloadGoal fnaFilesGoal,
			AdditionalFastasGoal additionalGoal,
			ObjectGoal<AccessionMap, GSProject> accessionMapGoal, FillSizeGoal fillSizeGoal,
			ObjectGoal<MurmurCGATBloomFilter, GSProject> bloomFilterGoal, Goal<GSProject>... deps) {
		super(project, name, project.getOutputFile(name, FileType.SER), Goal.append(deps, categoriesGoal,
				taxNodesGoal, fnaFilesGoal, accessionMapGoal, fillSizeGoal, bloomFilterGoal));
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
				getProject().getConfig().getKMerFastBloomFpp(), null, false);
		store.initSize(bloomFilterGoal.get().getEntries());

		try {
			Set<TaxIdNode> taxNodes = taxNodesGoal.get();
			MyFastaReader fastaReader = new MyFastaReader(getProject().getConfig().getMaxReadSizeBytes(),
					taxNodesGoal.get(), accessionMapGoal.get(), store);

			for (File fnaFile : fnaFilesGoal.getFiles()) {
				RefSeqCategory cat = fnaFilesGoal.getCategoryForFile(fnaFile);
				if (categoriesGoal.get()[0].contains(cat)) {
					fastaReader.readFasta(fnaFile);
				}
			}
			for (File additionalFasta : additionalGoal.getFiles()) {
				TaxIdNode additionalNode = additionalGoal.getTaxNodeForFile(additionalFasta);
				if (taxNodes.contains(additionalNode)) {
					fastaReader.ignoreAccessionMap(additionalNode);
					fastaReader.readFasta(additionalFasta);
				}
			}
			if (getLogger().isWarnEnabled()) {
				getLogger().warn("Not stored kmers: " + fastaReader.tooManyCounter);
			}

			store.optimize();

			KMerStoreWrapper wrapper = new KMerStoreWrapper((KMerSortedArray<String>) store);
			wrapper.save(storeFile);
			if (getLogger().isInfoEnabled()) {
				getLogger().info("File saved " + storeFile);
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
				KMerSortedArray<String> store) {
			super(bufferSize, taxNodes, accessionMap, store.getK());
			this.store = store;
		}

		@Override
		protected void handleStore() {
			if (store.isFull()) {
				tooManyCounter++;
			} else {
				store.put(byteRingBuffer, node.getTaxId(), false);
			}
		}
	}
}