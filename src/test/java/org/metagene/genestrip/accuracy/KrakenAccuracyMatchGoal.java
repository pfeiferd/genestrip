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
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.accuracy.AccuracyMatchGoal.AccuracyCounts;
import org.metagene.genestrip.goals.kraken.KrakenResCountGoal;
import org.metagene.genestrip.io.StreamProvider;
import org.metagene.genestrip.io.StreamingFileResource;
import org.metagene.genestrip.io.StreamingResource;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.tax.TaxTree;

public class KrakenAccuracyMatchGoal extends KrakenResCountGoal {
	private final ObjectGoal<TaxTree, GSProject> taxTreeGoal;
	private final AccuracyCounts accuracyCounts;

	@SafeVarargs
	public KrakenAccuracyMatchGoal(GSProject project, ObjectGoal<Map<String, List<StreamingResource>>, GSProject> fastqMapGoal,
			ObjectGoal<TaxTree, GSProject> taxTreeGoal, Goal<GSProject>... deps) {
		super(project, fastqMapGoal, null, append(deps, taxTreeGoal));
		this.taxTreeGoal = taxTreeGoal;
		accuracyCounts = new AccuracyCounts();
	}

	@Override
	protected void makeFile(File file) throws IOException {
		accuracyCounts.clear();
		List<File> files = new ArrayList<File>();
		for (StreamingResource s : fileToFastqs.get(file)) {
			files.add(((StreamingFileResource) s).getFile());
		}
		computeStats(files);
		try (PrintStream out = new PrintStream(StreamProvider.getOutputStreamForFile(file))) {
			accuracyCounts.printCounts(out);
		}
	}

	@Override
	protected void afterReadMatch(String krakenTaxid, byte[] readDescriptor) {
		accuracyCounts.updateCounts(krakenTaxid, readDescriptor, taxTreeGoal.get());
	}
}
