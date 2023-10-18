package org.metagene.genestrip.accuracy;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;

import org.metagene.genestrip.GSConfig;
import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.goals.MatchGoal;
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

public class AccuracyMatchGoal extends MatchGoal {
	private int taxIdCorrectCount;
	private int genusCorrectCount;
	private int totalCount;

	@SafeVarargs
	public AccuracyMatchGoal(GSProject project, String name, File fastq, ObjectGoal<TaxTree, GSProject> taxTreeGoal,
			ObjectGoal<KMerStoreWrapper, GSProject> storeGoal, Goal<GSProject>... deps) {
		super(project, name, fastq, taxTreeGoal, storeGoal, false, deps);
	}

	@Override
	protected FastqKMerMatcher createMatcher(KMerStoreWrapper wrapper, GSConfig config, TaxTree taxTree) {
		return new FastqKMerMatcher(wrapper.getKmerStore(), config.getMaxReadSizeBytes(), config.getThreadQueueSize(),
				config.getThreads(), config.getMaxKMerResCounts(), taxTree, config.getMaxReadTaxErrorCount()) {
			@Override
			protected boolean matchRead(MyReadEntry entry, boolean reverse) {
				totalCount++;
				boolean res = super.matchRead(entry, reverse);
				int colonIndex = ByteArrayUtil.indexOf(entry.readDescriptor, 1, entry.readDescriptorSize, ':');
				String correctTaxId = new String(entry.readDescriptor, 1, colonIndex);
				if (correctTaxId.equals(entry.readTaxId)) {
					taxIdCorrectCount++;
					genusCorrectCount++;
				}
				else if (entry.readTaxId != null && entry.readTaxId != INVALID_TAX) {
					TaxIdNode correctGenusTaxNode = getGenusTaxNode(correctTaxId);
					if (correctGenusTaxNode != null) {
						if (correctGenusTaxNode == getGenusTaxNode(entry.readTaxId)) {
							genusCorrectCount++;
						}
					}
				}
				return res;
			}
			
			protected TaxIdNode getGenusTaxNode(String taxid) {
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
		};
	}
	
	@Override
	protected void writeOutputFile(File file, MatchingResult result, KMerStoreWrapper wrapper) throws IOException {
		PrintStream out = new PrintStream(StreamProvider.getOutputStreamForFile(file));
		out.println("total; taxid correct; genus correct");
		out.print(totalCount);
		out.print(';');
		out.print(genusCorrectCount);
		out.print(';');
		out.print(taxIdCorrectCount);
		out.println(';');
		out.close();
	}
}
