package org.metagene.genestrip.goals;

import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.zip.GZIPInputStream;

import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.kraken.KrakenResultFastqMergeListener;
import org.metagene.genestrip.kraken.KrakenResultFastqMerger;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.tax.AssemblySummaryReader.FTPEntryQuality;
import org.metagene.genestrip.tax.AssemblySummaryReader.FTPEntryWithQuality;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;
import org.metagene.genestrip.trie.FastaTrieCleaner;
import org.metagene.genestrip.trie.KMerTrie;
import org.metagene.genestrip.util.CGATRingBuffer;

public class TrieFromKrakenResGoal extends ObjectGoal<KMerTrie<String>, GSProject> {
	private final ObjectGoal<Set<TaxIdNode>, GSProject> taxNodesGoal;
	private final ObjectGoal<Map<TaxIdNode, List<FTPEntryWithQuality>>, GSProject> fastaFilesGoal;

	@SafeVarargs
	public TrieFromKrakenResGoal(GSProject project, ObjectGoal<Set<TaxIdNode>, GSProject> taxNodesGoal,
			ObjectGoal<Map<TaxIdNode, List<FTPEntryWithQuality>>, GSProject> fastaFilesGoal,
			Goal<GSProject>... dependencies) {
		super(project, "triefromkrakenres", dependencies);
		this.taxNodesGoal = taxNodesGoal;
		this.fastaFilesGoal = fastaFilesGoal;
	}

	@Override
	public void makeThis() {
		try {
			KMerTrie<String> trie = new KMerTrie<String>(getProject().getkMserSize());
			Set<TaxIdNode> nodes = taxNodesGoal.get();
			byte[] kmer = new byte[trie.getLen()];

			KrakenResultFastqMerger merger = new KrakenResultFastqMerger(
					getProject().getConfig().getMaxReadSizeBytes());

			FileInputStream fStream = new FileInputStream(getProject().getFastqFile());
			GZIPInputStream gStream = new GZIPInputStream(fStream, 4096);

			merger.process(new BufferedInputStream(new FileInputStream(getProject().getKrakenOutFile())), gStream,
					KrakenResultFastqMergeListener.createFilterByTaxIdNodes(nodes,
							new KrakenResultFastqMergeListener() {

								@Override
								public void newTaxIdForRead(long lineCount, byte[] readDescriptor, byte[] read,
										byte[] readProbs, String krakenTaxid, int bps, String kmerTaxid, int hitLength,
										byte[] output) {
									for (int j = 0; j < hitLength; j++) {
										for (int i = 0; i < kmer.length; i++) {
											kmer[i] = read[post + j + i];
										}
										trie.put(kmer, 0, kmerTaxid);
									}
								}
							}));

			FastaTrieCleaner<String> fastaTrieCleaner = new FastaTrieCleaner<String>();

			FTPEntryQuality minQuality = getProject().getConfig().getFastaQuality();
			CGATRingBuffer ringBuffer = new CGATRingBuffer(trie.getLen());
			byte[] buffer = new byte[4096 * 8];

			for (TaxIdNode node : nodes) {
				List<FTPEntryWithQuality> entryList = fastaFilesGoal.get().get(node);

				for (FTPEntryWithQuality entry : entryList) {
					if (minQuality == null || !entry.getQuality().below(minQuality)) {

						File file = new File(getProject().getFastasDir(), entry.getFileName());
						if (file.exists()) {
							fastaTrieCleaner.fillIndex(trie, node.getTaxId(), file, buffer, ringBuffer);
						}
					}
				}
			}
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}
}
