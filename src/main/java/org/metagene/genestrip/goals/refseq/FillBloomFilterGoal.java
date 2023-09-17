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
import org.metagene.genestrip.util.ByteArrayUtil;

public class FillBloomFilterGoal extends ObjectGoal<MurmurCGATBloomFilter, GSProject> {
	private final Collection<RefSeqCategory> includedCategories;
	private final ObjectGoal<Set<TaxIdNode>, GSProject> taxNodesGoal;
	private final RefSeqFnaFilesDownloadGoal fnaFilesGoal;
	private final ObjectGoal<Map<String, TaxIdNode>, GSProject> accessionTrieGoal;
	private final ObjectGoal<Long, GSProject> sizeGoal;

	@SafeVarargs
	public FillBloomFilterGoal(GSProject project, String name, ObjectGoal<Set<TaxIdNode>, GSProject> taxNodesGoal,
			RefSeqFnaFilesDownloadGoal fnaFilesGoal, ObjectGoal<Map<String, TaxIdNode>, GSProject> accessionTrieGoal,
			FillSizeGoal sizeGoal, Goal<GSProject>... deps) {
		super(project, name, ArraysUtil.append(deps, taxNodesGoal, fnaFilesGoal, accessionTrieGoal, sizeGoal));
		this.includedCategories = sizeGoal.getIncludedCategories();
		this.taxNodesGoal = taxNodesGoal;
		this.fnaFilesGoal = fnaFilesGoal;
		this.accessionTrieGoal = accessionTrieGoal;
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
		private Map<String, TaxIdNode> accessionTrie;
		private Set<TaxIdNode> taxNodes;

		public MyFastaReader(AbstractKMerBloomIndex bloomIndex, int bufferSize) {
			super(bloomIndex, bufferSize);
			inCountRegion = false;
			accessionTrie = accessionTrieGoal.get();
			taxNodes = taxNodesGoal.get();
		}

		@Override
		protected void infoLine() throws IOException {
			super.infoLine();
			if (taxNodes.isEmpty()) {
				inCountRegion = true;
			} else {
				inCountRegion = false;
				int pos = ByteArrayUtil.indexOf(target, 0, size, ' ');
				if (pos >= 0) {
					String accession = new String(target, 1, pos - 1);
					TaxIdNode node = accessionTrie.get(accession);
					if (node != null) {
						inCountRegion = taxNodes.contains(node);
					}
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