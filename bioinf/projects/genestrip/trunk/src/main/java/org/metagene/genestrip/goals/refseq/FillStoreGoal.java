package org.metagene.genestrip.goals.refseq;

import java.io.File;
import java.io.IOException;
import java.util.Set;

import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.GSProject.FileType;
import org.metagene.genestrip.bloom.MurmurCGATBloomFilter;
import org.metagene.genestrip.make.FileListGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.store.KMerSortedArray;
import org.metagene.genestrip.store.KMerStoreWrapper;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;
import org.metagene.genestrip.util.ArraysUtil;

public class FillStoreGoal extends FileListGoal<GSProject> {
	private final ObjectGoal<Set<RefSeqCategory>[], GSProject> categoriesGoal;
	private final ObjectGoal<Set<TaxIdNode>, GSProject> taxNodesGoal;
	private final RefSeqFnaFilesDownloadGoal fnaFilesGoal;
	private final ObjectGoal<AccessionMap, GSProject> accessionMapGoal;
	private final ObjectGoal<MurmurCGATBloomFilter, GSProject> bloomFilterGoal;

	@SafeVarargs
	public FillStoreGoal(GSProject project, String name, ObjectGoal<Set<RefSeqCategory>[], GSProject> categoriesGoal, ObjectGoal<Set<TaxIdNode>, GSProject> taxNodesGoal,
			RefSeqFnaFilesDownloadGoal fnaFilesGoal, ObjectGoal<AccessionMap, GSProject> accessionMapGoal,
			FillSizeGoal fillSizeGoal, ObjectGoal<MurmurCGATBloomFilter, GSProject> bloomFilterGoal,
			Goal<GSProject>... deps) {
		super(project, name, project.getOutputFile(name, FileType.SER),
				ArraysUtil.append(deps, categoriesGoal, taxNodesGoal, fnaFilesGoal, accessionMapGoal, fillSizeGoal, bloomFilterGoal));
		this.categoriesGoal = categoriesGoal;
		this.taxNodesGoal = taxNodesGoal;
		this.fnaFilesGoal = fnaFilesGoal;
		this.accessionMapGoal = accessionMapGoal;
		this.bloomFilterGoal = bloomFilterGoal;
	}

	@Override
	public void makeFile(File storeFile) {
		KMerSortedArray<String> store = new KMerSortedArray<String>(getProject().getConfig().getKMerSize(), 0.000000001,
				null, false, false);
		store.initSize(bloomFilterGoal.get().getEntries());

		try {
			MyFastaReader fastaReader = new MyFastaReader(getProject().getConfig().getMaxReadSizeBytes(),
					taxNodesGoal.get(), accessionMapGoal.get(), store);

			for (File fnaFile : fnaFilesGoal.getFiles()) {
				RefSeqCategory cat = fnaFilesGoal.getCategoryForFile(fnaFile);
				if (categoriesGoal.get()[0].contains(cat)) {
					fastaReader.readFasta(fnaFile);
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