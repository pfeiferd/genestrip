package org.metagene.genestrip.goals.refseq;

import java.io.File;
import java.io.IOException;
import java.util.Set;

import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.bloom.MurmurCGATBloomFilter;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;
import org.metagene.genestrip.util.ArraysUtil;

public class FillBloomFilterGoal extends ObjectGoal<MurmurCGATBloomFilter, GSProject> {
	private final ObjectGoal<Set<RefSeqCategory>[], GSProject> categoriesGoal;
	private final ObjectGoal<Set<TaxIdNode>, GSProject> taxNodesGoal;
	private final RefSeqFnaFilesDownloadGoal fnaFilesGoal;
	private final ObjectGoal<AccessionMap, GSProject> accessionMapGoal;
	private final ObjectGoal<Long, GSProject> sizeGoal;

	@SafeVarargs
	public FillBloomFilterGoal(GSProject project, String name,
			ObjectGoal<Set<RefSeqCategory>[], GSProject> categoriesGoal,
			ObjectGoal<Set<TaxIdNode>, GSProject> taxNodesGoal, RefSeqFnaFilesDownloadGoal fnaFilesGoal,
			ObjectGoal<AccessionMap, GSProject> accessionMapGoal, FillSizeGoal sizeGoal, Goal<GSProject>... deps) {
		super(project, name, ArraysUtil.append(deps, taxNodesGoal, fnaFilesGoal, accessionMapGoal, sizeGoal));
		this.categoriesGoal = categoriesGoal;
		this.taxNodesGoal = taxNodesGoal;
		this.fnaFilesGoal = fnaFilesGoal;
		this.accessionMapGoal = accessionMapGoal;
		this.sizeGoal = sizeGoal;
	}

	@Override
	public void makeThis() {
		try {
			MurmurCGATBloomFilter filter = new MurmurCGATBloomFilter(getProject().getConfig().getKMerSize(),
					getProject().getConfig().getKMerFastBloomFpp());
			filter.ensureExpectedSize(sizeGoal.get(), false);

			MyFastaReader fastaReader = new MyFastaReader(getProject().getConfig().getMaxReadSizeBytes(),
					taxNodesGoal.get(), accessionMapGoal.get(), filter);

			for (File fnaFile : fnaFilesGoal.getFiles()) {
				RefSeqCategory cat = fnaFilesGoal.getCategoryForFile(fnaFile);
				if (categoriesGoal.get()[0].contains(cat)) {
					fastaReader.readFasta(fnaFile);
				}
			}

			set(filter);
			if (getLogger().isInfoEnabled()) {
				getLogger().info("Bloom filter entries: " + filter.getEntries());
			}
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}

	protected static class MyFastaReader extends AbstractStoreFastaReader {
		private final MurmurCGATBloomFilter filter;

		public MyFastaReader(int bufferSize, Set<TaxIdNode> taxNodes, AccessionMap accessionMap,
				MurmurCGATBloomFilter filter) {
			super(bufferSize, taxNodes, accessionMap, filter.getK());
			this.filter = filter;
		}

		@Override
		protected void handleStore() {
			if (!filter.contains(byteRingBuffer, false)) {
				filter.put(byteRingBuffer);
			}
		}
	}
}