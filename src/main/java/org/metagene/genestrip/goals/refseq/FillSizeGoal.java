package org.metagene.genestrip.goals.refseq;

import java.io.File;
import java.io.IOException;
import java.util.Set;

import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;
import org.metagene.genestrip.util.ArraysUtil;

public class FillSizeGoal extends ObjectGoal<Long, GSProject> {
	private final ObjectGoal<Set<RefSeqCategory>[], GSProject> categoriesGoal;
	private final ObjectGoal<Set<TaxIdNode>, GSProject> taxNodesGoal;
	private final RefSeqFnaFilesDownloadGoal fnaFilesGoal;
	private final ObjectGoal<AccessionMap, GSProject> accessionMapGoal;

	@SafeVarargs
	public FillSizeGoal(GSProject project, String name, ObjectGoal<Set<RefSeqCategory>[], GSProject> categoriesGoal,
			ObjectGoal<Set<TaxIdNode>, GSProject> taxNodesGoal, RefSeqFnaFilesDownloadGoal fnaFilesGoal,
			ObjectGoal<AccessionMap, GSProject> accessionMapGoal, Goal<GSProject>... deps) {
		super(project, name, ArraysUtil.append(deps, categoriesGoal, taxNodesGoal, fnaFilesGoal, accessionMapGoal));
		this.categoriesGoal = categoriesGoal;
		this.taxNodesGoal = taxNodesGoal;
		this.fnaFilesGoal = fnaFilesGoal;
		this.accessionMapGoal = accessionMapGoal;
	}

	@Override
	public void makeThis() {
		try {
			MyFastaReader fastaReader = new MyFastaReader(getProject().getConfig().getMaxReadSizeBytes(),
					taxNodesGoal.get(), accessionMapGoal.get(), getProject().getConfig().getKMerSize());

			for (File fnaFile : fnaFilesGoal.getFiles()) {
				RefSeqCategory cat = fnaFilesGoal.getCategoryForFile(fnaFile);
				if (categoriesGoal.get()[0].contains(cat)) {
					fastaReader.readFasta(fnaFile);
				}
			}
			if (getLogger().isInfoEnabled()) {
				getLogger().info("Store size determined: " + fastaReader.getCounter());
			}

			set(fastaReader.getCounter());
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}

	protected static class MyFastaReader extends AbstractRefSeqFastaReader {
		private long counter;
		private final int k;

		public MyFastaReader(int bufferSize, Set<TaxIdNode> taxNodes, AccessionMap accessionMap, int k) {
			super(bufferSize, taxNodes, accessionMap);
			this.k = k;
			counter = 0;
		}

		@Override
		protected void infoLine() throws IOException {
			if (includeRegion) {
				counter -= k - 1;
			}
			super.infoLine();
		}

		@Override
		protected void dataLine() throws IOException {
			if (includeRegion) {
				counter += size - 1;
			}
		}

		public long getCounter() {
			return counter;
		}

		@Override
		protected void done() throws IOException {
			if (includeRegion) {
				counter -= k - 1;
			}
		}
	}
}