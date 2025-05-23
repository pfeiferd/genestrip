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
package org.metagene.genestrip.goals.kraken;

import org.metagene.genestrip.GSConfigKey;
import org.metagene.genestrip.GSGoalKey;
import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.GSProject.FileType;
import org.metagene.genestrip.goals.MultiFileGoal;
import org.metagene.genestrip.io.StreamProvider;
import org.metagene.genestrip.io.StreamingFileResource;
import org.metagene.genestrip.io.StreamingResource;
import org.metagene.genestrip.io.StreamingResourceStream;
import org.metagene.genestrip.kraken.KrakenExecutor;
import org.metagene.genestrip.kraken.KrakenResultListener;
import org.metagene.genestrip.kraken.KrakenResultProcessor;
import org.metagene.genestrip.make.FileListGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.match.MatchingResult;
import org.metagene.genestrip.match.ResultReporter;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;
import org.metagene.genestrip.util.DigitTrie;

import java.io.*;
import java.util.*;

public class KrakenResFileGoal extends FileListGoal<GSProject> {
	private final ObjectGoal<Map<String, StreamingResourceStream>, GSProject> fastqMapGoal;
	private final ObjectGoal<Map<String, List<KrakenResCountGoal.KrakenResStats>>, GSProject> countGoal;
	private final Map<File, String> fileToKeyMap;

	@SafeVarargs
	public KrakenResFileGoal(GSProject project,
                             ObjectGoal<Map<String, StreamingResourceStream>, GSProject> fastqMapGoal,
                             ObjectGoal<Set<TaxIdNode>, GSProject> taxNodesGoal,
							 ObjectGoal<Map<String, List<KrakenResCountGoal.KrakenResStats>>, GSProject> countGoal,
							 Goal<GSProject>... deps) {
		super(project, GSGoalKey.KRAKENRES, (List<File>) null, Goal.append(deps, taxNodesGoal, countGoal, fastqMapGoal));
		this.fastqMapGoal = fastqMapGoal;
		this.countGoal = countGoal;
		fileToKeyMap = new HashMap<>();
	}

	@Override
	protected void provideFiles() {
		for (String key : fastqMapGoal.get().keySet()) {
			File resFile = getProject().getOutputFile(getKey().getName(), key, null, FileType.KRAKEN_OUT_RES, false);
			addFile(resFile);
			fileToKeyMap.put(resFile, key);
		}
	}

	@Override
	protected void makeFile(File file) throws IOException {
		List<KrakenResCountGoal.KrakenResStats> list = countGoal.get().get(fileToKeyMap.get(file));
		try (PrintStream out = new PrintStream(StreamProvider.getOutputStreamForFile(file))) {
			print(list, out);
		}
	}

	protected void print(List<KrakenResCountGoal.KrakenResStats> list, PrintStream out) {
		out.println("taxid;reads;kmers;kmers in matching reads");
		for (KrakenResCountGoal.KrakenResStats stats : list) {
			out.print(stats.getTaxid());
			out.print(';');
			out.print(stats.getReads());
			out.print(';');
			out.print(stats.getKmers());
			out.print(';');
			out.print(stats.getKmersInMatchingReads());
			out.println(';');
		}
	}
}
