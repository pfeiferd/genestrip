package org.metagene.genestrip.goals.refseq;

import java.io.File;
import java.io.IOException;
import java.util.Collection;

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
	private final Collection<RefSeqCategory> categories;
	private final RefSeqFnaFilesDownloadGoal fnaFilesGoal;
	private final ObjectGoal<AccessionMap, GSProject> accessionTrieGoal;
	private final FillStoreGoal includeStoreGoal;
	private final ObjectGoal<TaxTree, GSProject> taxTreeGoal;

	@SafeVarargs
	public UpdateStoreGoal(GSProject project, String name, ObjectGoal<TaxTree, GSProject> taxTreeGoal,
			RefSeqFnaFilesDownloadGoal fnaFilesGoal, ObjectGoal<AccessionMap, GSProject> accessionTrieGoal,
			FillStoreGoal includeStoreGoal, Goal<GSProject>... deps) {
		super(project, name, project.getOutputFile(name, FileType.SER),
				ArraysUtil.append(deps, taxTreeGoal, fnaFilesGoal, accessionTrieGoal, includeStoreGoal));
		this.categories = fnaFilesGoal.getCategories();
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
				if (categories.contains(cat)) {
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