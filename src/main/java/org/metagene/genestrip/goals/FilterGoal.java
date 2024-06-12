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
import java.util.List;
import java.util.Map;

import org.metagene.genestrip.ExecutionContext;
import org.metagene.genestrip.GSConfigKey;
import org.metagene.genestrip.GSGoalKey;
import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.GSProject.FileType;
import org.metagene.genestrip.bloom.FastqBloomFilter;
import org.metagene.genestrip.io.StreamingResource;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;

public class FilterGoal extends MultiFileGoal {
	private final BloomIndexedGoal indexedGoal;
	private final ExecutionContext executorServiceBundle;

	@SafeVarargs
	public FilterGoal(GSProject project, 
			ObjectGoal<Map<String, List<StreamingResource>>, GSProject> fastqMapGoal, BloomIndexedGoal indexedGoal,
			ExecutionContext executorServiceBundle, Goal<GSProject>... deps) {
		super(project, GSGoalKey.FILTER, fastqMapGoal, append(deps, indexedGoal));
		this.executorServiceBundle = executorServiceBundle;
		this.indexedGoal = indexedGoal;
	}

	@Override
	protected FileType getFileType() {
		return FileType.FASTQ_RES;
	}

	@Override
	protected void makeFile(File file) {
		FastqBloomFilter f = null;
		try {
			GSProject project = getProject();
			List<StreamingResource> resources = fileToFastqs.get(file);
			File dumpFile = booleanConfigValue(GSConfigKey.WRITED_DUMPED_FASTQ)
					? project.getOutputFile("dumped", file.getName(), FileType.FASTQ_RES)
					: null;
			f = new FastqBloomFilter(indexedGoal.get(), intConfigValue(GSConfigKey.MIN_POS_COUNT_FILTER),
					doubleConfigValue(GSConfigKey.POS_RATIO_FILTER),
					intConfigValue(GSConfigKey.INITIAL_READ_SIZE_BYTES), intConfigValue(GSConfigKey.THREAD_QUEUE_SIZE),
					executorServiceBundle);
			f.runFilter(resources, file, dumpFile);
		} catch (IOException e) {
			throw new RuntimeException(e);
		} finally {
			if (f != null) {
				f.dump();
			}
		}
	}
	
	@Override
	protected boolean isUseGZip() {
		return true;
	}
}
