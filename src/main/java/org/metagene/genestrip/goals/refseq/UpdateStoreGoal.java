package org.metagene.genestrip.goals.refseq;

import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;
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
import org.metagene.genestrip.util.CGATRingBuffer;

public class UpdateStoreGoal extends FileListGoal<GSProject> {
	private final Collection<RefSeqCategory> categories;
	private final ObjectGoal<Set<TaxIdNode>, GSProject> taxNodesGoal;
	private final RefSeqFnaFilesDownloadGoal fnaFilesGoal;
	private final AccessionCollectionGoal accessionCollectionGoal;
	private final IncludeStoreGoal includeStoreGoal;
	private final ObjectGoal<TaxTree, GSProject> taxTreeGoal;

	@SafeVarargs
	public UpdateStoreGoal(GSProject project, String name, Collection<RefSeqCategory> categories,
			ObjectGoal<TaxTree, GSProject> taxTreeGoal, ObjectGoal<Set<TaxIdNode>, GSProject> taxNodesGoal,
			RefSeqFnaFilesDownloadGoal fnaFilesGoal, AccessionCollectionGoal accessionCollectionGoal,
			IncludeStoreGoal includeStoreGoal, Goal<GSProject>... deps) {
		super(project, name, project.getOutputFile(name, FileType.SER),
				ArraysUtil.append(deps, taxTreeGoal,taxNodesGoal, fnaFilesGoal, accessionCollectionGoal, includeStoreGoal));
		this.categories = categories;
		this.taxTreeGoal = taxTreeGoal;
		this.taxNodesGoal = taxNodesGoal;
		this.fnaFilesGoal = fnaFilesGoal;
		this.accessionCollectionGoal = accessionCollectionGoal;
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
		private final Map<String, TaxIdNode> accessionMap;
		private final KMerSortedArray<String> store;
		private final CGATRingBuffer byteRingBuffer;

		private final UpdateValueProvider<String> provider;

		public MyFastaReader(int kmerSize, int bufferSize, KMerSortedArray<String> store) {
			super(bufferSize);
			taxTree = taxTreeGoal.get();
			accessionMap = accessionCollectionGoal.get();
			byteRingBuffer = new CGATRingBuffer(kmerSize);
			this.store = store;
			provider = new UpdateValueProvider<String>() {
				@Override
				public String getUpdateValue(String oldValue) {
					TaxIdNode oldNode = taxTree.getNodeByTaxId(oldValue);
					
					TaxIdNode newNode = taxTree.getLeastCommonAncestor(oldNode, node);
					if (newNode == oldNode) {
						return oldValue;
					}
					return newNode.getTaxId();
				}
			};
		}

		@Override
		protected void infoLine() throws IOException {
			int i = 0;
			for (; i < size; i++) {
				if (target[i] == ' ' || target[i] == '\n') {
					break;
				}
			}
			String accession = new String(target, 1, i);
			node = accessionMap.get(accession);
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