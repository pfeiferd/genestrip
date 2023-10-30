package org.metagene.genestrip.accuracy;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;

import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.accuracy.AccuracyMatchGoal.AccuracyCounts;
import org.metagene.genestrip.goals.kraken.KrakenResCountGoal;
import org.metagene.genestrip.io.StreamProvider;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.tax.TaxTree;

public class KrakenAccuracyMatchGoal extends KrakenResCountGoal {
	private final ObjectGoal<TaxTree, GSProject> taxTreeGoal;
	private final AccuracyCounts accuracyCounts;

	@SafeVarargs
	public KrakenAccuracyMatchGoal(GSProject project, String name, File csvFile,
			ObjectGoal<TaxTree, GSProject> taxTreeGoal, Goal<GSProject>... deps) {
		super(project, name, true, csvFile, null, append(deps, taxTreeGoal));
		this.taxTreeGoal = taxTreeGoal;
		accuracyCounts = new AccuracyCounts();
	}

	@Override
	protected void makeFile(File file) throws IOException {
		accuracyCounts.clear();
		computeStats(fileToFastqs.get(file));
		PrintStream out = new PrintStream(StreamProvider.getOutputStreamForFile(file));
		accuracyCounts.printCounts(out);
		out.close();
	}

	@Override
	protected void afterReadMatch(String krakenTaxid, byte[] readDescriptor) {
		accuracyCounts.updateCounts(krakenTaxid, readDescriptor, taxTreeGoal.get());
	}
}
