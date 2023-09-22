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
import org.metagene.genestrip.make.FileListGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.store.KMerSortedArray;
import org.metagene.genestrip.store.KMerSortedArray.UpdateValueProvider;
import org.metagene.genestrip.store.KMerStoreWrapper;
import org.metagene.genestrip.tax.TaxTree;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;
import org.metagene.genestrip.util.ArraysUtil;

public class UpdateStoreGoal extends FileListGoal<GSProject> {
	private final ObjectGoal<Set<RefSeqCategory>[], GSProject> categoriesGoal;
	private final RefSeqFnaFilesDownloadGoal fnaFilesGoal;
	private final ObjectGoal<AccessionMap, GSProject> accessionTrieGoal;
	private final FillStoreGoal includeStoreGoal;
	private final ObjectGoal<TaxTree, GSProject> taxTreeGoal;

	@SafeVarargs
	public UpdateStoreGoal(GSProject project, String name, ObjectGoal<Set<RefSeqCategory>[], GSProject> categoriesGoal,
			ObjectGoal<TaxTree, GSProject> taxTreeGoal, RefSeqFnaFilesDownloadGoal fnaFilesGoal,
			ObjectGoal<AccessionMap, GSProject> accessionTrieGoal, FillStoreGoal includeStoreGoal,
			Goal<GSProject>... deps) {
		super(project, name, project.getOutputFile(name, FileType.SER), ArraysUtil.append(deps, categoriesGoal,
				taxTreeGoal, fnaFilesGoal, accessionTrieGoal, includeStoreGoal));
		this.categoriesGoal = categoriesGoal;
		this.taxTreeGoal = taxTreeGoal;
		this.fnaFilesGoal = fnaFilesGoal;
		this.accessionTrieGoal = accessionTrieGoal;
		this.includeStoreGoal = includeStoreGoal;
	}

	@Override
	public void makeFile(File storeFile) {
		try {
			KMerStoreWrapper wrapper = KMerStoreWrapper.load(includeStoreGoal.getFile());
			KMerSortedArray<String> store = wrapper.getKmerStore();

			MyFastaReader fastaReader = new MyFastaReader(getProject().getConfig().getMaxReadSizeBytes(),
					taxTreeGoal.get(), accessionTrieGoal.get(), store);

			for (File fnaFile : fnaFilesGoal.getFiles()) {
				RefSeqCategory cat = fnaFilesGoal.getCategoryForFile(fnaFile);
				if (categoriesGoal.get()[0].contains(cat)) {
					fastaReader.readFasta(fnaFile);
				}
			}

			System.out.println("** Hit rate: " + (taxTreeGoal.get().hits / (double) taxTreeGoal.get().total));

			KMerStoreWrapper wrapper2 = new KMerStoreWrapper((KMerSortedArray<String>) store);
			wrapper2.save(storeFile);
			if (getLogger().isInfoEnabled()) {
				getLogger().info("File saved " + storeFile);
			}
		} catch (IOException e) {
			throw new RuntimeException(e);
		} catch (ClassNotFoundException e) {
			throw new RuntimeException(e);
		}
	}

	protected class MyFastaReader extends AbstractStoreFastaReader {
		private final KMerSortedArray<String> store;

		private final UpdateValueProvider<String> provider;

		public MyFastaReader(int bufferSize, TaxTree taxTree, AccessionMap accessionMap,
				KMerSortedArray<String> store) {
			super(bufferSize, null, accessionMap, store.getK());
			this.store = store;
			provider = new UpdateValueProvider<String>() {
				@Override
				public String getUpdateValue(String oldValue) {
					TaxIdNode oldNode = taxTree.getNodeByTaxId(oldValue);
					TaxIdNode newNode = taxTree.getLeastCommonAncestor(oldNode, node);
					if (newNode == null || newNode == oldNode) {
						return oldValue;
					}
					return newNode.getTaxId();
				}
			};
			includeRegion = true;
		}

		@Override
		protected void infoLine() throws IOException {
			byteRingBuffer.reset();
			updateNodeFromInfoLine();
		}

		@Override
		protected void handleStore() {
			store.update(byteRingBuffer, provider, false);
		}
	}
}