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
import org.metagene.genestrip.match.FastqKMerMatcher;
import org.metagene.genestrip.match.MatchingResult;
import org.metagene.genestrip.store.KMerStoreWrapper;
import org.metagene.genestrip.tax.TaxTree;
import org.metagene.genestrip.tax.TaxTree.Rank;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;
import org.metagene.genestrip.util.ByteArrayUtil;

public class AccuracyMatchGoal extends MultiMatchGoal {
	private Map<String, Integer> noTaxIdErrorPerTaxid;
	private int taxIdCorrectCount;
	private int genusCorrectCount;
	private int genusIncorrectCount;
	private int noTaxIdCount;
	private int totalCount;
	private long startMillis;
	private boolean accuracy;

	@SafeVarargs
	public AccuracyMatchGoal(GSProject project, String name, File csvFile, ObjectGoal<TaxTree, GSProject> taxTreeGoal,
			ObjectGoal<KMerStoreWrapper, GSProject> storeGoal, Goal<GSProject>... deps) {
		super(project, name, csvFile, taxTreeGoal, storeGoal, false, deps);
	}

	@Override
	protected void makeFile(File file) {
		taxIdCorrectCount = 0;
		genusCorrectCount = 0;
		genusIncorrectCount = 0;
		noTaxIdCount = 0;
		totalCount = 0;
		noTaxIdErrorPerTaxid = new HashMap<String, Integer>();
		startMillis = System.currentTimeMillis();
		accuracy = file.getName().contains("accuracy");
		super.makeFile(file);
	}

	@Override
	protected FastqKMerMatcher createMatcher(KMerStoreWrapper wrapper, GSConfig config, TaxTree taxTree) {
		return new FastqKMerMatcher(wrapper.getKmerStore(), config.getMaxReadSizeBytes(), config.getThreadQueueSize(),
				config.getThreads(), config.getMaxKMerResCounts(), taxTree, config.getMaxReadTaxErrorCount()) {
			@Override
			protected void afterMatch(MyReadEntry entry, boolean found) throws IOException {
//				ByteArrayUtil.print(entry.readDescriptor, System.out);
				super.afterMatch(entry, found);
				totalCount++;
				if (accuracy) {
					int colonIndex = ByteArrayUtil.indexOf(entry.readDescriptor, 1, entry.readDescriptorSize, ':');
					String correctTaxId = new String(entry.readDescriptor, 1, colonIndex - 1);
					if (correctTaxId.equals(entry.readTaxId)) {
						taxIdCorrectCount++;
						genusCorrectCount++;
					} else if (entry.readTaxId != null && entry.readTaxId != INVALID_TAX) {
						TaxIdNode correctGenusTaxNode = taxTree.getRankedNode(correctTaxId, Rank.GENUS);
						if (correctGenusTaxNode != null) {
							if (correctGenusTaxNode == taxTree.getRankedNode(entry.readTaxId, Rank.GENUS)) {
								genusCorrectCount++;
							} else {
								genusIncorrectCount++;
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
		};
	}

	@Override
	protected void writeOutputFile(File file, MatchingResult result, KMerStoreWrapper wrapper) throws IOException {
		PrintStream out = new PrintStream(StreamProvider.getOutputStreamForFile(file));
		if (accuracy) {
			out.println("total; taxid correct; genus correct; genus incorrect; no taxid;");
			out.print(totalCount);
			out.print(';');
			out.print(taxIdCorrectCount);
			out.print(';');
			out.print(genusCorrectCount);
			out.print(';');
			out.print(genusIncorrectCount);
			out.print(';');
			out.print(noTaxIdCount);
			out.println(';');
			System.out.println(noTaxIdErrorPerTaxid);
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
}
