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
import java.util.List;
import java.util.Map;

import org.metagene.genestrip.ExecutionContext;
import org.metagene.genestrip.GSConfigKey;
import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.goals.MatchGoal;
import org.metagene.genestrip.io.StreamProvider;
import org.metagene.genestrip.io.StreamingResource;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.GoalKey;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.match.FastqKMerMatcher;
import org.metagene.genestrip.match.MatchingResult;
import org.metagene.genestrip.store.Database;
import org.metagene.genestrip.store.KMerSortedArray;
import org.metagene.genestrip.tax.Rank;
import org.metagene.genestrip.tax.SmallTaxTree;
import org.metagene.genestrip.tax.SmallTaxTree.SmallTaxIdNode;
import org.metagene.genestrip.tax.TaxTree;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;
import org.metagene.genestrip.util.ByteArrayUtil;

public class AccuracyMatchGoal extends MatchGoal {
	private final ObjectGoal<TaxTree, GSProject> taxTreeGoal;
	private final AccuracyCounts accuracyCounts;
	private long startMillis;
	private boolean timing;
	private int totalCount;

	@SafeVarargs
	public AccuracyMatchGoal(GSProject project, GoalKey goalKey,
			ObjectGoal<Map<String, List<StreamingResource>>, GSProject> fastqMapGoal,
			ObjectGoal<Database, GSProject> storeGoal, ObjectGoal<TaxTree, GSProject> taxTreeGoal,
			ExecutionContext bundle, Goal<GSProject>... deps) {
		super(project, goalKey, fastqMapGoal, storeGoal, bundle, append(deps, taxTreeGoal));
		this.taxTreeGoal = taxTreeGoal;
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

	@Override
	protected FastqKMerMatcher createMatcher(KMerSortedArray<SmallTaxIdNode> store, SmallTaxTree taxTree,
			ExecutionContext bundle) {
		return new FastqKMerMatcher(store, intConfigValue(GSConfigKey.INITIAL_READ_SIZE_BYTES),
				intConfigValue(GSConfigKey.THREAD_QUEUE_SIZE), bundle, intConfigValue(GSConfigKey.MAX_KMER_RES_COUNTS),
				booleanConfigValue(GSConfigKey.CLASSIFY_READS) ? taxTree : null,
				intConfigValue(GSConfigKey.MAX_CLASSIFICATION_PATHS),
				intConfigValue(GSConfigKey.MAX_READ_TAX_ERROR_COUNT)) {
			@Override
			protected void afterMatch(MyReadEntry entry, boolean found) throws IOException {
				super.afterMatch(entry, found);
				totalCount++;
				if (!timing) {
					String taxid = entry.classNode == null ? null : entry.classNode.getTaxId();
					accuracyCounts.updateCounts(taxid, entry.readDescriptor, taxTreeGoal.get());
				}
			}
		};
	}

	@Override
	protected void writeOutputFile(File file, MatchingResult result, Database wrapper) throws IOException {
		try (PrintStream out = new PrintStream(StreamProvider.getOutputStreamForFile(file))) {
			if (!timing) {
				accuracyCounts.printCounts(out);
				accuracyCounts.printNoTaxidErrors(System.out);
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
		}
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

		public void printNoTaxidErrors(PrintStream out) {
			out.println("taxid; error;");
			for (String tax : noTaxIdErrorPerTaxid.keySet()) {
				out.print(tax);
				out.print(";");
				out.print(noTaxIdErrorPerTaxid.get(tax));
				out.println(";");
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
