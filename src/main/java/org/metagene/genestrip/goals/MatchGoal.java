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

import org.metagene.genestrip.GSConfig;
import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.GSProject.FileType;
import org.metagene.genestrip.make.FileListGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.store.FastqKMerMatcher;
import org.metagene.genestrip.store.FastqKMerMatcher.Result;
import org.metagene.genestrip.store.KMerStoreWrapper;
import org.metagene.genestrip.store.ResultReporter;
import org.metagene.genestrip.util.ArraysUtil;
import org.metagene.genestrip.util.StreamProvider;

public class MatchGoal extends FileListGoal<GSProject> {
	private final File fastq;
	private final KMerStoreFileGoal storeGoal;
	private final boolean writedFiltered;

	@SafeVarargs
	public MatchGoal(GSProject project, String name, File fastq, KMerStoreFileGoal storeGoal, boolean writeFiltered,
			Goal<GSProject>... deps) {
		super(project, name, project.getOutputFile(name, fastq, FileType.CSV, false),
				ArraysUtil.append(deps, storeGoal));
		this.fastq = fastq;
		this.storeGoal = storeGoal;
		this.writedFiltered = writeFiltered;
	}

	@Override
	protected void makeFile(File file) {
		FastqKMerMatcher c = null;
		try {
			File filteredFile = null;
			File krakenOutStyleFile = null;
			if (writedFiltered) {
				filteredFile = getProject().getOutputFile(getName(), fastq, FileType.FASTQ_RES, true);
				krakenOutStyleFile = getProject().getOutputFile(getName(), fastq, FileType.KRAKEN_OUT_RES, false);
			}

			KMerStoreWrapper wrapper = KMerStoreWrapper.load(storeGoal.getFile());

			GSConfig config = getProject().getConfig();
			c = new FastqKMerMatcher(wrapper.getKmerStore(), config.getMaxReadSizeBytes(),
					config.getThreadQueueSize(), config.getThreads(), config.isCountUniqueKmers());
			Result res = c.runClassifier(fastq, filteredFile, krakenOutStyleFile);
			c.dump();

			PrintStream out = new PrintStream(StreamProvider.getOutputStreamForFile(file));
			new ResultReporter(wrapper.getTaxids()).print(res, out);
			out.close();
		} catch (IOException | ClassNotFoundException e) {
			throw new RuntimeException(e);
		} finally {
			if (c != null) {
				c.dump();
			}
		}
	}	
}
