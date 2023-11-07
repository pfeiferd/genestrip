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
package org.metagene.genestrip.accuracy;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.HashMap;
import java.util.Map;

import org.metagene.genestrip.GSConfig;
import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.goals.MultiMatchGoal;
import org.metagene.genestrip.io.StreamProvider;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.match.FastqKMerMatcher2;
import org.metagene.genestrip.match.MatchingResult;
import org.metagene.genestrip.store.KMerSortedArray;
import org.metagene.genestrip.store.KMerSortedArray.ValueConverter;
import org.metagene.genestrip.store.KMerStoreWrapper;
import org.metagene.genestrip.tax.TaxTree;
import org.metagene.genestrip.tax.TaxTree.Rank;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;
import org.metagene.genestrip.util.ByteArrayUtil;

public class AccuracyMatchGoal extends MultiMatchGoal {
	private final AccuracyCounts accuracyCounts;
	private long startMillis;
	private boolean timing;
	private int totalCount;

	@SafeVarargs
	public AccuracyMatchGoal(GSProject project, String name, File csvFile, ObjectGoal<TaxTree, GSProject> taxTreeGoal,
			ObjectGoal<KMerStoreWrapper, GSProject> storeGoal, Goal<GSProject>... deps) {
		super(project, name, true, csvFile, taxTreeGoal, storeGoal, false, deps);
		accuracyCounts = new AccuracyCounts();
	}

	@Override
	protected void makeFile(File file) {
		accuracyCounts.clear();
		startMillis = System.currentTimeMillis();
		totalCount = 0;
		timing = file.getName().contains("timing");
		super.makeFile(file);
	}

	protected FastqKMerMatcher2 createMatcher(KMerStoreWrapper wrapper, TaxTree taxTree) {
		GSConfig config = getProject().getConfig();

		KMerSortedArray<TaxIdNode> store = new KMerSortedArray<TaxIdNode>(wrapper.getKmerStore(),
				new ValueConverter<String, TaxIdNode>() {
					@Override
					public TaxIdNode convertValue(String value) {
						return taxTree.getNodeByTaxId(value);
					}
				});

		return new FastqKMerMatcher2(store, config.getMaxReadSizeBytes(), config.getThreadQueueSize(),
				config.getThreads(), config.getMaxKMerResCounts(), getProject().isClassifyReads() ? taxTree : null,
				config.getMaxClassificationPaths(), getProject().getMaxReadTaxErrorCount());
	}
	
	protected FastqKMerMatcher2 createMatcher(KMerSortedArray<TaxIdNode> store, TaxTree taxTree) {
		GSConfig config = getProject().getConfig();

		return new FastqKMerMatcher2(store, config.getMaxReadSizeBytes(), config.getThreadQueueSize(),
				config.getThreads(), config.getMaxKMerResCounts(), getProject().isClassifyReads() ? taxTree : null,
				config.getMaxClassificationPaths(), getProject().getMaxReadTaxErrorCount()) {
			@Override
			protected void afterMatch(MyReadEntry entry, boolean found) throws IOException {
				super.afterMatch(entry, found);
				totalCount++;
				if (!timing) {
					String taxid = entry.classNode == null ? null : entry.classNode.getTaxId();
					accuracyCounts.updateCounts(taxid, entry.readDescriptor, taxTree);
				}
			}
		};

	}

	@Override
	protected void writeOutputFile(File file, MatchingResult result, KMerStoreWrapper wrapper) throws IOException {
		PrintStream out = new PrintStream(StreamProvider.getOutputStreamForFile(file));
		if (!timing) {
			accuracyCounts.printCounts(out);
		} else {
			long millis = System.currentTimeMillis() - startMillis;
			out.println("total; elapsed millis; reads per min.;");
			out.print(totalCount);
			out.print(';');
			out.print(millis);
			out.print(';');
			out.print(totalCount * 1000 * 60 / millis);
			out.println(';');
		}
		out.close();
	}

	public static class AccuracyCounts {
		public final Map<String, Integer> noTaxIdErrorPerTaxid;
		public final Map<Rank, Integer> correctRankOnPath;

		public int taxIdCorrectCount;
		public int genusCorrectCount;
		public int genusIncorrectCount;
		public int noTaxIdCount;
		public long totalCount;

		public AccuracyCounts() {
			noTaxIdErrorPerTaxid = new HashMap<String, Integer>();
			correctRankOnPath = new HashMap<Rank, Integer>();
		}

		public void clear() {
			taxIdCorrectCount = 0;
			genusCorrectCount = 0;
			genusIncorrectCount = 0;
			noTaxIdCount = 0;
			totalCount = 0;

			noTaxIdErrorPerTaxid.clear();
			correctRankOnPath.clear();
		}

		public void updateCounts(String readTaxId, byte[] readDescriptor, TaxTree taxTree) {
			totalCount++;
			int colonIndex = ByteArrayUtil.indexOf(readDescriptor, 1, readDescriptor.length, ':');
			if (colonIndex == -1) {
				ByteArrayUtil.print(readDescriptor, System.out);
			} else {
				String correctTaxId = new String(readDescriptor, 1, colonIndex - 1);
				if (readTaxId != null) {
					if (correctTaxId.equals(readTaxId)) {
						taxIdCorrectCount++;
						genusCorrectCount++;
					} else {
						TaxIdNode correctGenusTaxNode = taxTree.getRankedNode(correctTaxId, Rank.GENUS);
						if (correctGenusTaxNode != null) {
							if (correctGenusTaxNode == taxTree.getRankedNode(readTaxId, Rank.GENUS)) {
								genusCorrectCount++;
							} else {
								genusIncorrectCount++;
							}
						}
					}

					for (Rank rank : Rank.values()) {
						TaxIdNode correctRankTaxNode = taxTree.getRankedNode(correctTaxId, rank);
						if (correctRankTaxNode != null) {
							TaxIdNode node = taxTree.getNodeByTaxId(readTaxId);
							if (correctRankTaxNode.equals(node)) {
								Integer c = correctRankOnPath.get(rank);
								if (c == null) {
									c = 0;
								}
								correctRankOnPath.put(rank, c + 1);
							}
						}
					}
				} else {
					Integer e = noTaxIdErrorPerTaxid.get(correctTaxId);
					if (e == null) {
						e = 0;
					}
					noTaxIdErrorPerTaxid.put(correctTaxId, e + 1);
					noTaxIdCount++;
				}
			}
		}

		public void printCounts(PrintStream out) {
			out.print("total; taxid correct; genus correct; genus incorrect; no taxid;");
			for (Rank rank : Rank.values()) {
				if (correctRankOnPath.get(rank) != null) {
					out.print(rank.getName() + " correct on path;");
				}
			}
			out.println();
			out.print(totalCount);
			out.print(';');
			out.print(taxIdCorrectCount);
			out.print(';');
			out.print(genusCorrectCount);
			out.print(';');
			out.print(genusIncorrectCount);
			out.print(';');
			out.print(noTaxIdCount);
			out.print(';');
			for (Rank rank : Rank.values()) {
				Integer c = correctRankOnPath.get(rank);
				if (c != null) {
					out.print(c);
					out.print(';');
				}
			}
			out.println();
		}
	}
}
