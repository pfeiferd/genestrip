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
import java.nio.charset.StandardCharsets;

import org.metagene.genestrip.GSGoalKey;
import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.GSProject.GSFileType;
import org.metagene.genestrip.io.StreamProvider;
import org.metagene.genestrip.io.StreamingResourceStream;
import org.metagene.genestrip.make.FileListGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;

import java.io.*;
import java.util.*;

/**
 * Writes the per-taxid Kraken statistics for each key to a CSV output file.
 *
 * @param <P> the project type
 */
public class KrakenResFileGoal<P extends GSProject> extends FileListGoal<P> {
	private final ObjectGoal<Map<String, StreamingResourceStream>, P> fastqMapGoal;
	private final ObjectGoal<Map<String, List<KrakenResCountGoal.KrakenResStats>>, P> countGoal;
	private final Map<File, String> fileToKeyMap;

	/**
	 * Creates the goal that writes the Kraken result files.
	 *
	 * @param project      the project type
	 * @param fastqMapGoal the goal providing the FASTQ streams keyed by name
	 * @param taxNodesGoal the goal providing the tax id nodes of interest
	 * @param countGoal    the goal providing the per-key Kraken result statistics
	 * @param deps         additional goals this goal depends on
	 */
	@SafeVarargs
	public KrakenResFileGoal(P project,
                             ObjectGoal<Map<String, StreamingResourceStream>, P> fastqMapGoal,
                             ObjectGoal<Set<TaxIdNode>, P> taxNodesGoal,
							 ObjectGoal<Map<String, List<KrakenResCountGoal.KrakenResStats>>, P> countGoal,
							 Goal<P>... deps) {
		super(project, GSGoalKey.KRAKENRES, (List<File>) null, Goal.append(deps, taxNodesGoal, countGoal, fastqMapGoal));
		this.fastqMapGoal = fastqMapGoal;
		this.countGoal = countGoal;
		fileToKeyMap = new HashMap<>();
	}

	@Override
	protected void provideFiles() {
		for (String key : fastqMapGoal.get().keySet()) {
			File resFile = getProject().getOutputFile(getKey().getName(), key, null, GSFileType.KRAKEN_OUT_RES, false);
			addFile(resFile);
			fileToKeyMap.put(resFile, key);
		}
	}

	@Override
	protected void makeFile(File file) throws IOException {
		List<KrakenResCountGoal.KrakenResStats> list = countGoal.get().get(fileToKeyMap.get(file));
		try (PrintStream out = new PrintStream(StreamProvider.getOutputStreamForFile(file), false, StandardCharsets.UTF_8)) {
			print(list, out);
		}
	}

	/**
	 * Writes the given per-taxid statistics as CSV (taxid, reads, kmers, kmers in matching reads) to {@code out}.
	 *
	 * @param list the per-taxid statistics to write
	 * @param out  the stream to write the CSV rows to
	 */
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
