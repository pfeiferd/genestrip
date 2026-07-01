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
import java.util.Map;

import org.metagene.genestrip.ExecutionContext;
import org.metagene.genestrip.GSConfigKey;
import org.metagene.genestrip.GSGoalKey;
import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.GSProject.GSFileType;
import org.metagene.genestrip.bloom.FastqBloomFilter;
import org.metagene.genestrip.io.StreamingResourceStream;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;

/**
 * Filters each input FASTQ file against the loaded bloom-filter index, writing the reads that pass to a per-key
 * FASTQ output file (and optionally the rejected reads to a dump file).
 *
 * @param <P> the project type
 */
public class FilterGoal<P extends GSProject> extends MultiFileGoal<P> {
	private final LoadIndexGoal<P> indexedGoal;
	private final ExecutionContext executorServiceBundle;

	/**
	 * Creates the goal, wiring the input FASTQ map, the loaded index and the execution context.
	 *
	 * @param project the project type
	 * @param fastqMapGoal the goal supplying the per-key input FASTQ resources
	 * @param indexedGoal the goal supplying the loaded bloom-filter index
	 * @param executorServiceBundle the execution context providing threading and shared services
	 * @param deps the additional goals this goal depends on
	 */
	@SafeVarargs
	public FilterGoal(P project,
			ObjectGoal<Map<String, StreamingResourceStream>, P> fastqMapGoal, LoadIndexGoal<P> indexedGoal,
			ExecutionContext executorServiceBundle, Goal<P>... deps) {
		super(project, GSGoalKey.FILTER, fastqMapGoal, append(deps, indexedGoal));
		this.executorServiceBundle = executorServiceBundle;
		this.indexedGoal = indexedGoal;
	}

	@Override
	protected GSFileType getFileType() {
		return GSFileType.FASTQ_RES;
	}

	@Override
	protected boolean isUseGZip() {
		return booleanConfigValue(GSConfigKey.GZIP_FASTQ_OUTPUT);
	}

	@Override
	protected void makeFile(File file) throws IOException {
		FastqBloomFilter f = null;
		try {
			P project = getProject();
			StreamingResourceStream resources = fileToFastqs.get(file);
			File dumpFile = booleanConfigValue(GSConfigKey.WRITE_DUMPED_FASTQ)
					? project.getOutputFile("dumped", null, file.getName(), GSFileType.FASTQ_RES, isUseGZip())
					: null;
			f = new FastqBloomFilter(intConfigValue(GSConfigKey.KMER_SIZE), indexedGoal.get(), intConfigValue(GSConfigKey.MIN_POS_COUNT_FILTER),
					doubleConfigValue(GSConfigKey.POS_RATIO_FILTER),
					intConfigValue(GSConfigKey.INITIAL_READ_SIZE_BYTES), intConfigValue(GSConfigKey.THREAD_QUEUE_SIZE),
					executorServiceBundle, booleanConfigValue(GSConfigKey.WITH_PROBS)) {
				@Override
				protected boolean isProgressBar() {
					return booleanConfigValue(GSConfigKey.PROGRESS_BAR);
				}

				@Override
				protected String getProgressBarTaskName() {
					return getKey().getName();
				}
			};
			f.runFilter(resources, file, dumpFile);
		} finally {
			if (f != null) {
				f.dump();
			}
		}
	}
}
