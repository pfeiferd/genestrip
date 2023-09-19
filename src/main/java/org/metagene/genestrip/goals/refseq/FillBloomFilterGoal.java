package org.metagene.genestrip.goals.refseq;

import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.Set;

import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.bloom.MurmurCGATBloomFilter;
import org.metagene.genestrip.fasta.AbstractFastaReader;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;
import org.metagene.genestrip.util.ArraysUtil;
import org.metagene.genestrip.util.ByteArrayUtil;
import org.metagene.genestrip.util.CGAT;
import org.metagene.genestrip.util.CGATRingBuffer;

public class FillBloomFilterGoal extends ObjectGoal<MurmurCGATBloomFilter, GSProject> {
	private final Collection<RefSeqCategory> includedCategories;
	private final ObjectGoal<Set<TaxIdNode>, GSProject> taxNodesGoal;
	private final RefSeqFnaFilesDownloadGoal fnaFilesGoal;
	private final ObjectGoal<AccessionMap, GSProject> accessionTrieGoal;
	private final ObjectGoal<Long, GSProject> sizeGoal;

	@SafeVarargs
	public FillBloomFilterGoal(GSProject project, String name, ObjectGoal<Set<TaxIdNode>, GSProject> taxNodesGoal,
			RefSeqFnaFilesDownloadGoal fnaFilesGoal, ObjectGoal<AccessionMap, GSProject> accessionTrieGoal,
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
			MurmurCGATBloomFilter filter = new MurmurCGATBloomFilter(getProject().getConfig().getKMerSize(),
					0.000000001);
			filter.ensureExpectedSize(sizeGoal.get(), false);

			MyFastaReader fastaReader = new MyFastaReader(getProject().getConfig().getMaxReadSizeBytes(), filter);

			for (File fnaFile : fnaFilesGoal.getFiles()) {
				RefSeqCategory cat = fnaFilesGoal.getCategoryForFile(fnaFile);
				if (includedCategories.contains(cat)) {
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

	protected class MyFastaReader extends AbstractFastaReader {
		private boolean inStoreRegion;
		private TaxIdNode node;
		private final AccessionMap accessionTrie;
		private final Set<TaxIdNode> taxNodes;
		private final MurmurCGATBloomFilter filter;
		private final CGATRingBuffer byteRingBuffer;

		public MyFastaReader(int bufferSize, MurmurCGATBloomFilter filter) {
			super(bufferSize);
			inStoreRegion = false;
			accessionTrie = accessionTrieGoal.get();
			taxNodes = taxNodesGoal.get();
			byteRingBuffer = new CGATRingBuffer(filter.getK());
			this.filter = filter;
		}

		@Override
		protected void infoLine() throws IOException {
			inStoreRegion = false;
			int pos = ByteArrayUtil.indexOf(target, 0, size, ' ');
			if (pos >= 0) {
				node = accessionTrie.get(target, 1, pos);
				if (node != null) {
					inStoreRegion = taxNodes.isEmpty() || taxNodes.contains(node);
				}
			}
		}

		@Override
		protected void dataLine() {
			if (inStoreRegion) {
				for (int i = 0; i < size - 1; i++) {
					byteRingBuffer.put(CGAT.cgatToUpperCase(target[i]));
					if (byteRingBuffer.isFilled() && byteRingBuffer.isCGAT()) {
						if (!filter.contains(byteRingBuffer, false)) {
							filter.put(byteRingBuffer);
						}
					}
				}
			}
		}
	}
}