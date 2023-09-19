package org.metagene.genestrip.goals.refseq;

import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.Set;

import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.GSProject.FileType;
import org.metagene.genestrip.bloom.MurmurCGATBloomFilter;
import org.metagene.genestrip.fasta.AbstractFastaReader;
import org.metagene.genestrip.make.FileListGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.store.KMerSortedArray;
import org.metagene.genestrip.store.KMerStoreWrapper;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;
import org.metagene.genestrip.util.ArraysUtil;
import org.metagene.genestrip.util.ByteArrayUtil;
import org.metagene.genestrip.util.CGAT;
import org.metagene.genestrip.util.CGATRingBuffer;

public class FillStoreGoal extends FileListGoal<GSProject> {
	private final Collection<RefSeqCategory> includeCategories;
	private final ObjectGoal<Set<TaxIdNode>, GSProject> taxNodesGoal;
	private final RefSeqFnaFilesDownloadGoal fnaFilesGoal;
	private final ObjectGoal<AccessionMap, GSProject> accessionTrieGoal;
	private final ObjectGoal<MurmurCGATBloomFilter, GSProject> bloomFilterGoal;

	@SafeVarargs
	public FillStoreGoal(GSProject project, String name, ObjectGoal<Set<TaxIdNode>, GSProject> taxNodesGoal,
			RefSeqFnaFilesDownloadGoal fnaFilesGoal, ObjectGoal<AccessionMap, GSProject> accessionTrieGoal,
			FillSizeGoal fillSizeGoal, ObjectGoal<MurmurCGATBloomFilter, GSProject> bloomFilterGoal,
			Goal<GSProject>... deps) {
		super(project, name, project.getOutputFile(name, FileType.SER),
				ArraysUtil.append(deps, taxNodesGoal, fnaFilesGoal, accessionTrieGoal, fillSizeGoal, bloomFilterGoal));
		this.includeCategories = fillSizeGoal.getIncludedCategories();
		this.taxNodesGoal = taxNodesGoal;
		this.fnaFilesGoal = fnaFilesGoal;
		this.accessionTrieGoal = accessionTrieGoal;
		this.bloomFilterGoal = bloomFilterGoal;
	}

	@Override
	public void makeFile(File storeFile) {
		KMerSortedArray<String> store = new KMerSortedArray<String>(getProject().getConfig().getKMerSize(), 0.000000001,
				null, false, false, bloomFilterGoal.get());
		store.initSize(bloomFilterGoal.get().getEntries());

		try {
			MyFastaReader fastaReader = new MyFastaReader(getProject().getConfig().getMaxReadSizeBytes(), store);

			for (File fnaFile : fnaFilesGoal.getFiles()) {
				RefSeqCategory cat = fnaFilesGoal.getCategoryForFile(fnaFile);
				if (includeCategories.contains(cat)) {
					fastaReader.readFasta(fnaFile);
				}
			}
			if (getLogger().isInfoEnabled()) {
				getLogger().info("Not stored kmers: " + fastaReader.tooManyCounter);
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

	protected class MyFastaReader extends AbstractFastaReader {
		private boolean inStoreRegion;
		private TaxIdNode node;
		private final AccessionMap accessionTrie;
		private final Set<TaxIdNode> taxNodes;
		private final KMerSortedArray<String> store;
		private final CGATRingBuffer byteRingBuffer;

		private long tooManyCounter;

		public MyFastaReader(int bufferSize, KMerSortedArray<String> store) {
			super(bufferSize);
			inStoreRegion = false;
			accessionTrie = accessionTrieGoal.get();
			taxNodes = taxNodesGoal.get();
			byteRingBuffer = new CGATRingBuffer(store.getK());
			this.store = store;
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
		protected void dataLine() throws IOException {
			if (inStoreRegion) {
				for (int i = 0; i < size - 1; i++) {
					byteRingBuffer.put(CGAT.cgatToUpperCase(target[i]));
					if (byteRingBuffer.isFilled() && byteRingBuffer.isCGAT()) {
						if (store.isFull()) {
							tooManyCounter++;
						} else {
							store.put(byteRingBuffer, node.getTaxId(), false);
						}
					}
				}
			}
		}
	}
}