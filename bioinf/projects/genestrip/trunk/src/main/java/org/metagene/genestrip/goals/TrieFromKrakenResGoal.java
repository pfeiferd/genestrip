/*
 * 
 * “Commons Clause” License Condition v1.0
 * 
 * The Software is provided to you by the Licensor under the License, 
 * as defined below, subject to the following condition.
 * 
 * Without limiting other conditions in the License, the grant of rights under the License 
 * will not include, and the License does not grant to you, the right to Sell the Software.
 * 
 * For purposes of the foregoing, “Sell” means practicing any or all of the rights granted 
 * to you under the License to provide to third parties, for a fee or other consideration 
 * (including without limitation fees for hosting or consulting/ support services related to 
 * the Software), a product or service whose value derives, entirely or substantially, from the 
 * functionality of the Software. Any license notice or attribution required by the License 
 * must also include this Commons Clause License Condition notice.
 * 
 * Software: genestrip
 * 
 * License: Apache 2.0
 * 
 * Licensor: Daniel Pfeifer (daniel.pfeifer@progotec.de)
 * 
 */
package org.metagene.genestrip.goals;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.Serializable;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.goals.TrieFromKrakenResGoal.TaxidWithCount;
import org.metagene.genestrip.kraken.KrakenResultFastqMergeListener;
import org.metagene.genestrip.kraken.KrakenResultFastqMerger;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.tax.AssemblySummaryReader.FTPEntryQuality;
import org.metagene.genestrip.tax.AssemblySummaryReader.FTPEntryWithQuality;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;
import org.metagene.genestrip.trie.FastaTrieCleaner;
import org.metagene.genestrip.trie.KMerTrie;
import org.metagene.genestrip.trie.KMerTrie.KMerTrieVisitor;
import org.metagene.genestrip.util.ByteArrayUtil;
import org.metagene.genestrip.util.CGATRingBuffer;
import org.metagene.genestrip.util.StreamProvider;

public class TrieFromKrakenResGoal extends ObjectGoal<KMerTrie<TaxidWithCount>, GSProject> {
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
			KMerTrie<TaxidWithCount> trie = new KMerTrie<TaxidWithCount>(getProject().getkMserSize());
			Set<TaxIdNode> nodes = taxNodesGoal.get();
			byte[] kmer = new byte[trie.getLen()];

			KrakenResultFastqMerger merger = new KrakenResultFastqMerger(
					getProject().getConfig().getMaxReadSizeBytes());

			KrakenResultFastqMergeListener printListener = KrakenResultFastqMergeListener
					.createPrintListener(System.out, null);

			InputStream stream1 = StreamProvider.getInputStreamForFile(getProject().getKrakenOutFile());
			InputStream stream2 = StreamProvider.getInputStreamForFile(getProject().getFastqFile());			
			merger.process(stream1, stream2, KrakenResultFastqMergeListener
							.createFilterByTaxIdNodes(nodes, new KrakenResultFastqMergeListener() {
								private long lastLineCount = -1;

								@Override
								public void newTaxIdForRead(long lineCount, byte[] readDescriptor, byte[] read,
										byte[] readProbs, String krakenTaxid, int bps, int pos, String kmerTaxid,
										int hitLength, byte[] output) {
									if (lastLineCount != lineCount) {
										lastLineCount = lineCount;
										printListener.newTaxIdForRead(lineCount, readDescriptor, read, readProbs,
												krakenTaxid, bps, pos, kmerTaxid, hitLength, output);
									}
									for (int j = 0; j < hitLength; j++) {
										for (int i = 0; i < kmer.length; i++) {
											kmer[i] = read[pos + j + i];
										}
										System.out.println("Position:  " + pos);
										System.out.println(kmerTaxid);
										ByteArrayUtil.print(kmer, System.out);
										System.out.println();

										TaxidWithCount tc = trie.get(kmer, 0, false);
										if (tc == null) {
											tc = new TaxidWithCount(kmerTaxid);
											trie.put(kmer, 0, tc, false);
										}
										tc.inc();
									}
								}
							}));
			stream1.close();
			stream2.close();

			FTPEntryQuality minQuality = getProject().getConfig().getFastaQuality();
			String[] matchingTaxId = new String[1];
			CGATRingBuffer ringBuffer = new CGATRingBuffer(trie.getLen());
			
			FastaTrieCleaner<TaxidWithCount> fastaTrieCleaner = new FastaTrieCleaner<TaxidWithCount>(trie, ringBuffer, 4096) {
				@Override
				protected boolean isMatchingValue(TaxidWithCount value) {
					return matchingTaxId[0].equals(value.taxid);
				}
			};

			for (TaxIdNode node : nodes) {
				matchingTaxId[0] = node.getTaxId();

				List<FTPEntryWithQuality> entryList = fastaFilesGoal.get().get(node);

				if (entryList != null) {
					for (FTPEntryWithQuality entry : entryList) {
						if (minQuality == null || !entry.getQuality().below(minQuality)) {

							File file = new File(getProject().getFastasDir(), entry.getFileName());
							if (getLogger().isInfoEnabled()) {
								getLogger().info("Cleaning via file " + file);
							}
							if (file.exists()) {
								fastaTrieCleaner.readFasta(file);
							}
						}
					}
				}
			}

			trie.visit(new KMerTrieVisitor<TaxidWithCount>() {
				@Override
				public void nextValue(KMerTrie<TaxidWithCount> trie, byte[] kmer, TaxidWithCount value) {
					ByteArrayUtil.print(kmer, System.out);
					System.out.println();
				}
			}, false);

			set(trie);
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}

	@SuppressWarnings("serial")
	public static class TaxidWithCount implements Serializable {
		private int count;
		private final String taxid;

		public TaxidWithCount(String taxid) {
			this.taxid = taxid;
			count = 0;
		}

		public void inc() {
			count++;
		}

		public int getCount() {
			return count;
		}

		public String getTaxid() {
			return taxid;
		}
	}
}
