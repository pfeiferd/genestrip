package org.metagene.genestrip.goals.refseq;

import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;
import java.util.Set;

import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.bloom.AbstractKMerBloomIndex;
import org.metagene.genestrip.bloom.FastaIndexer;
import org.metagene.genestrip.bloom.MurmurCGATBloomFilter;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;
import org.metagene.genestrip.util.ArraysUtil;

public class IncludeBloomFilterGoal extends ObjectGoal<MurmurCGATBloomFilter, GSProject> {
	private final Collection<RefSeqCategory> includedCategories;
	private final ObjectGoal<Set<TaxIdNode>, GSProject> taxNodesGoal;
	private final RefSeqFnaFilesDownloadGoal fnaFilesGoal;
	private final ObjectGoal<Map<String, TaxIdNode>, GSProject> accessionCollectionGoal;
	private final ObjectGoal<Long, GSProject> sizeGoal;

	@SafeVarargs
	public IncludeBloomFilterGoal(GSProject project, String name, ObjectGoal<Set<TaxIdNode>, GSProject> taxNodesGoal,
			RefSeqFnaFilesDownloadGoal fnaFilesGoal,
			ObjectGoal<Map<String, TaxIdNode>, GSProject> accessionCollectionGoal, IncludeSizeGoal sizeGoal,
			Goal<GSProject>... deps) {
		super(project, name, ArraysUtil.append(deps, taxNodesGoal, fnaFilesGoal, accessionCollectionGoal, sizeGoal));
		this.includedCategories = sizeGoal.getIncludedCategories();
		this.taxNodesGoal = taxNodesGoal;
		this.fnaFilesGoal = fnaFilesGoal;
		this.accessionCollectionGoal = accessionCollectionGoal;
		this.sizeGoal = sizeGoal;
	}

	@Override
	public void makeThis() {
		try {
			AbstractKMerBloomIndex bloomIndex = new AbstractKMerBloomIndex("temp",
					getProject().getConfig().getKMerSize(), 0.000000001, null);
			bloomIndex.getFilter().ensureExpectedSize(sizeGoal.get(), false);

			MyFastaReader fastaReader = new MyFastaReader(bloomIndex, getProject().getConfig().getMaxReadSizeBytes());

			for (File fnaFile : fnaFilesGoal.getFiles()) {
				RefSeqCategory cat = fnaFilesGoal.getCategoryForFile(fnaFile);
				if (includedCategories.contains(cat)) {
					fastaReader.readFasta(fnaFile);
				}
			}

			set(bloomIndex.getFilter());
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}

	protected class MyFastaReader extends FastaIndexer {
		private boolean inCountRegion;
		private Map<String, TaxIdNode> accessionMap;
		private Set<TaxIdNode> taxNodes;

		public MyFastaReader(AbstractKMerBloomIndex bloomIndex, int bufferSize) {
			super(bloomIndex, bufferSize);
			inCountRegion = false;
			accessionMap = accessionCollectionGoal.get();
			taxNodes = taxNodesGoal.get();
		}

		@Override
		protected void infoLine() throws IOException {
			super.infoLine();
			if (taxNodes.isEmpty()) {
				inCountRegion = true;
			} else {
				int i = 0;
				for (; i < size; i++) {
					if (target[i] == ' ' || target[i] == '\n') {
						break;
					}
				}
				String accession = new String(target, 1, i);
				TaxIdNode node = accessionMap.get(accession);
				if (node != null) {
					inCountRegion = taxNodes.contains(node);
				} else {
					inCountRegion = false;
				}
			}
		}

		@Override
		protected void dataLine() {
			if (inCountRegion) {
				super.dataLine();
			}
		}
	}
}