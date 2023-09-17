package org.metagene.genestrip.goals.refseq;

import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.Collections;
import java.util.Map;
import java.util.Set;

import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.fasta.AbstractFastaReader;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;
import org.metagene.genestrip.util.ArraysUtil;
import org.metagene.genestrip.util.ByteArrayUtil;

public class FillSizeGoal extends ObjectGoal<Long, GSProject> {
	private final Collection<RefSeqCategory> includedCategories;
	private final ObjectGoal<Set<TaxIdNode>, GSProject> taxNodesGoal;
	private final RefSeqFnaFilesDownloadGoal fnaFilesGoal;
	private final ObjectGoal<Map<String, TaxIdNode>, GSProject> accessionCollectionGoal;

	private final int kmerSize;

	@SafeVarargs
	public FillSizeGoal(GSProject project, String name, Collection<RefSeqCategory> includeCategories,
			ObjectGoal<Set<TaxIdNode>, GSProject> taxNodesGoal, RefSeqFnaFilesDownloadGoal fnaFilesGoal,
			ObjectGoal<Map<String, TaxIdNode>, GSProject> accessionCollectionGoal, Goal<GSProject>... deps) {
		super(project, name, ArraysUtil.append(deps, taxNodesGoal, fnaFilesGoal, accessionCollectionGoal));
		this.includedCategories = Collections.unmodifiableCollection(includeCategories);
		this.taxNodesGoal = taxNodesGoal;
		this.fnaFilesGoal = fnaFilesGoal;
		this.accessionCollectionGoal = accessionCollectionGoal;
		kmerSize = project.getConfig().getKMerSize();
	}

	public Collection<RefSeqCategory> getIncludedCategories() {
		return includedCategories;
	}

	@Override
	public void makeThis() {
		try {
			MyFastaReader fastaReader = new MyFastaReader(getProject().getConfig().getMaxReadSizeBytes());

			for (File fnaFile : fnaFilesGoal.getFiles()) {
				RefSeqCategory cat = fnaFilesGoal.getCategoryForFile(fnaFile);
				if (includedCategories.contains(cat)) {
					fastaReader.readFasta(fnaFile);
				}
			}

			set(fastaReader.getCounter());
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}

	protected class MyFastaReader extends AbstractFastaReader {
		private long counter;
		private boolean inCountRegion;
		private Map<String, TaxIdNode> accessionTrie;
		private Set<TaxIdNode> taxNodes;

		public MyFastaReader(int bufferSize) {
			super(bufferSize);
			counter = 0;
			inCountRegion = false;
			accessionTrie = accessionCollectionGoal.get();
			taxNodes = taxNodesGoal.get();
		}

		@Override
		protected void infoLine() throws IOException {
			if (inCountRegion) {
				counter -= kmerSize - 1;
			}
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
		protected void dataLine() throws IOException {
			if (inCountRegion) {
				counter += size - 1;
			}
		}

		public long getCounter() {
			return counter;
		}

		@Override
		protected void done() throws IOException {
			counter -= kmerSize + 1;
		}
	}
}