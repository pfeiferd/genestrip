package org.metagene.genestrip.goals.refseq;

import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.Set;

import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.GSProject.FileType;
import org.metagene.genestrip.fasta.AbstractFastaReader;
import org.metagene.genestrip.make.FileListGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.store.KMerSortedArray;
import org.metagene.genestrip.store.KMerSortedArray.UpdateValueProvider;
import org.metagene.genestrip.store.KMerStoreWrapper;
import org.metagene.genestrip.tax.TaxTree;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;
import org.metagene.genestrip.util.ArraysUtil;
import org.metagene.genestrip.util.ByteArrayUtil;
import org.metagene.genestrip.util.CGATRingBuffer;

public class UpdateStoreGoal extends FileListGoal<GSProject> {
	private final Collection<RefSeqCategory> categories;
	private final ObjectGoal<Set<TaxIdNode>, GSProject> taxNodesGoal;
	private final RefSeqFnaFilesDownloadGoal fnaFilesGoal;
	private final ObjectGoal<AccessionMap, GSProject> accessionTrieGoal;
	private final FillStoreGoal includeStoreGoal;
	private final ObjectGoal<TaxTree, GSProject> taxTreeGoal;

	@SafeVarargs
	public UpdateStoreGoal(GSProject project, String name, ObjectGoal<TaxTree, GSProject> taxTreeGoal,
			ObjectGoal<Set<TaxIdNode>, GSProject> taxNodesGoal, RefSeqFnaFilesDownloadGoal fnaFilesGoal,
			ObjectGoal<AccessionMap, GSProject> accessionTrieGoal, FillStoreGoal includeStoreGoal,
			Goal<GSProject>... deps) {
		super(project, name, project.getOutputFile(name, FileType.SER), ArraysUtil.append(deps, taxTreeGoal,
				taxNodesGoal, fnaFilesGoal, accessionTrieGoal, includeStoreGoal));
		this.categories = fnaFilesGoal.getCategories();
		this.taxTreeGoal = taxTreeGoal;
		this.taxNodesGoal = taxNodesGoal;
		this.fnaFilesGoal = fnaFilesGoal;
		this.accessionTrieGoal = accessionTrieGoal;
		this.includeStoreGoal = includeStoreGoal;
	}

	@Override
	public void makeFile(File storeFile) {
		try {
			KMerStoreWrapper wrapper = KMerStoreWrapper.load(includeStoreGoal.getFile());
			KMerSortedArray<String> store = wrapper.getKmerStore();

			MyFastaReader fastaReader = new MyFastaReader(getProject().getConfig().getKMerSize(),
					getProject().getConfig().getMaxReadSizeBytes(), store);

			for (File fnaFile : fnaFilesGoal.getFiles()) {
				RefSeqCategory cat = fnaFilesGoal.getCategoryForFile(fnaFile);
				if (categories.contains(cat)) {
					fastaReader.readFasta(fnaFile);
				}
			}

			// TODO: How to best create the store stats?
			KMerStoreWrapper wrapper2 = new KMerStoreWrapper((KMerSortedArray<String>) store, taxNodesGoal.get(), null);
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

	protected class MyFastaReader extends AbstractFastaReader {
		private TaxIdNode node;

		private final TaxTree taxTree;
		private final AccessionMap accessionTrie;
		private final KMerSortedArray<String> store;
		private final CGATRingBuffer byteRingBuffer;

		private final UpdateValueProvider<String> provider;

		public MyFastaReader(int kmerSize, int bufferSize, KMerSortedArray<String> store) {
			super(bufferSize);
			taxTree = taxTreeGoal.get();
			accessionTrie = accessionTrieGoal.get();
			byteRingBuffer = new CGATRingBuffer(kmerSize);
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
		}

		@Override
		protected void infoLine() throws IOException {
			node = null;
			int pos = ByteArrayUtil.indexOf(target, 0, size, ' ');
			if (pos >= 0) {
				node = accessionTrie.get(target, 1, pos);
			}
		}

		@Override
		protected void dataLine() throws IOException {
			for (int i = 0; i < size - 1; i++) {
				byteRingBuffer.put(target[i]);
				if (byteRingBuffer.isFilled() && byteRingBuffer.isCGAT()) {
					store.update(byteRingBuffer, provider, false);
				}
			}
		}
	}
}