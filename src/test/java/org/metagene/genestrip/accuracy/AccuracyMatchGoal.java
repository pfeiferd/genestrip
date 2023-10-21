package org.metagene.genestrip.accuracy;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;

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
	private int taxIdCorrectCount;
	private int genusCorrectCount;
	private int genusIncorrectCount;
	private int noTaxIdCount;
	private int totalCount;

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
				int colonIndex = ByteArrayUtil.indexOf(entry.readDescriptor, 1, entry.readDescriptorSize, ':');
				String correctTaxId = new String(entry.readDescriptor, 1, colonIndex - 1);
				if (correctTaxId.equals(entry.readTaxId)) {
					taxIdCorrectCount++;
					genusCorrectCount++;
				}
				else if (entry.readTaxId != null && entry.readTaxId != INVALID_TAX) {
					TaxIdNode correctGenusTaxNode = getGenusTaxNode(taxTree, correctTaxId);
					if (correctGenusTaxNode != null) {
						if (correctGenusTaxNode == getGenusTaxNode(taxTree, entry.readTaxId)) {
							genusCorrectCount++;
						}
						else {
							genusIncorrectCount++;
						}
					}
				}
				else {
					noTaxIdCount++;
				}
			}
		};
	}
	
	public static TaxIdNode getGenusTaxNode(TaxTree taxTree, String taxid) {
		TaxIdNode node = taxTree.getNodeByTaxId(taxid);
		while (node != null) {
			if (node.getRank() == Rank.GENUS) {
				return node;
			}
			if (!node.getRank().isBelow(Rank.GENUS)) {
				return null;
			}
			node = node.getParent();
		}
		return null;
	}
	
	@Override
	protected void writeOutputFile(File file, MatchingResult result, KMerStoreWrapper wrapper) throws IOException {
		PrintStream out = new PrintStream(StreamProvider.getOutputStreamForFile(file));
		out.println("total; taxid correct; genus correct; genus inccorrect; no taxid;");
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
		out.close();
	}
}
