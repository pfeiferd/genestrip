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
import java.io.PrintStream;
import java.util.List;

import org.metagene.genestrip.GSConfig;
import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.GSProject.FileType;
import org.metagene.genestrip.make.FileListGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.store.FastqClassifier;
import org.metagene.genestrip.store.KMerStore;
import org.metagene.genestrip.store.FastqClassifier.StatsPerTaxid;
import org.metagene.genestrip.util.ArraysUtil;
import org.metagene.genestrip.util.StreamProvider;

public class ClassifyGoal extends FileListGoal<GSProject> {
	private final File fastq;
	private final KMerStoreFileGoal trieGoal;
	private final boolean writedFiltered;

	@SafeVarargs
	public ClassifyGoal(GSProject project, String name, File fastq, KMerStoreFileGoal trieGoal, boolean writeFiltered,
			Goal<GSProject>... deps) {
		super(project, name, project.getOutputFile(name, fastq, FileType.CSV, false),
				ArraysUtil.append(deps, trieGoal));
		this.fastq = fastq;
		this.trieGoal = trieGoal;
		this.writedFiltered = writeFiltered;
	}

	@Override
	protected void makeFile(File file) {
		FastqClassifier c = null;
		try {
			File filteredFile = null;
			File krakenOutStyleFile = null;
			if (writedFiltered) {
				filteredFile = getProject().getOutputFile(getName(), fastq, FileType.FASTQ_RES, true);
				krakenOutStyleFile = getProject().getOutputFile(getName(), fastq, FileType.KRAKEN_OUT_RES, false);
			}

			KMerStore<String> trie = KMerStore.load(trieGoal.getFile());

			GSConfig config = getProject().getConfig();
			c = new FastqClassifier(trie, config.getMaxReadSizeBytes(),
					config.getThreadQueueSize(), config.getThreads(), config.isCountUniqueKmers());
			List<StatsPerTaxid> res = c.runClassifier(fastq, filteredFile, krakenOutStyleFile);
			c.dump();

			PrintStream out = new PrintStream(StreamProvider.getOutputStreamForFile(file));
			print(res, out);
			out.close();
		} catch (IOException | ClassNotFoundException e) {
			throw new RuntimeException(e);
		} finally {
			if (c != null) {
				c.dump();
			}
		}
	}
	
	public static void print(List<StatsPerTaxid> allStats, PrintStream out) {
		out.println("taxid;reads;kmers;unique kmers;contigs;average contig length;max contig length;");
		for (StatsPerTaxid stats : allStats) {
			out.print(stats.getTaxid());
			out.print(';');
			out.print(stats.getReads());
			out.print(';');
			out.print(stats.getKmers());
			out.print(';');
			out.print(stats.getUniqueKmers());
			out.print(';');
			out.print(stats.getContigs());
			out.print(';');
			out.print(((double) stats.getKmers()) / stats.getContigs());
			out.print(';');
			out.print(stats.getMaxContigLen());
			out.println(';');
		}		
	}
}