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

import org.metagene.genestrip.GSConfig;
import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.GSProject.FileType;
import org.metagene.genestrip.bloom.FastqBloomFilter;
import org.metagene.genestrip.goals.refseq.BloomIndexedGoal;
import org.metagene.genestrip.make.FileListGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.util.ExecutorServiceBundle;

public class FilterGoal extends FileListGoal<GSProject> {
	private final File fastq;
	private final BloomIndexedGoal indexedGoal;
	private final ExecutorServiceBundle executorServiceBundle;

	@SafeVarargs
	public FilterGoal(GSProject project, String name, File fastq, BloomIndexedGoal indexedGoal,
			ExecutorServiceBundle executorServiceBundle, Goal<GSProject>... deps) {
		super(project, name, project.getOutputFile("filtered", fastq.getName(), FileType.FASTQ_RES), append(deps, indexedGoal));
		this.executorServiceBundle = executorServiceBundle;
		this.fastq = fastq;
		this.indexedGoal = indexedGoal;
	}

	@Override
	protected void makeFile(File file) {
		FastqBloomFilter f = null;
		try {
			GSProject project = getProject();
			GSConfig config = project.getConfig();

			File dumpFile = config.isWriteDumpedFastq() ? project.getOutputFile("dumped", fastq.getName(), FileType.FASTQ_RES)
					: null;
			f = new FastqBloomFilter(indexedGoal.get(), config.getMinPosCountFilter(), config.getPosRatioFilter(),
					config.getMaxReadSizeBytes(), config.getThreadQueueSize(), executorServiceBundle);
			f.runFilter(fastq, file, dumpFile);
		} catch (IOException e) {
			throw new RuntimeException(e);
		} finally {
			if (f != null) {
				f.dump();
			}
		}
	}
}
