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
import java.util.LinkedHashMap;
import java.util.Map;

import org.metagene.genestrip.GSGoalKey;
import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.GSProject.FileType;
import org.metagene.genestrip.io.StreamingFileResource;
import org.metagene.genestrip.io.StreamingResource;
import org.metagene.genestrip.io.StreamingResourceListStream;
import org.metagene.genestrip.io.StreamingResourceStream;
import org.metagene.genestrip.io.StreamingURLResource;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;

public class FastqMapTransformGoal extends ObjectGoal<Map<String, StreamingResourceStream>, GSProject> {
	private final ObjectGoal<Map<String, StreamingResourceStream>, GSProject> inputGoal;

	private final boolean fastqType;

	@SafeVarargs
	public FastqMapTransformGoal(GSProject project, boolean fastqType,
			ObjectGoal<Map<String, StreamingResourceStream>, GSProject> inputGoal, Goal<GSProject>... deps) {
		super(project, fastqType ? GSGoalKey.FASTQ_MAP_TRANSFORM : GSGoalKey.FASTA_MAP_TRANSFORM, append(deps, inputGoal));
		this.inputGoal = inputGoal;
		this.fastqType = fastqType;
	}

	@Override
	protected void doMakeThis() {
		if (!fastqType || getProject().isDownloadFastqs() || getProject().isDownloadFastqsToCommon()) {
			// Linked hash map preserve order of keys as entered.
			Map<String, StreamingResourceStream> map = new LinkedHashMap<String, StreamingResourceStream>();
			if (!fastqType || getProject().isDownloadFastqsToCommon()) {
				fillMap(map, fastqType ? getProject().getCommon().getFastqDir() : getProject().getCommon().getFastaDir());
			}
			else if (!fastqType || getProject().isDownloadFastqs()) {
				fillMap(map, fastqType ? getProject().getFastqDir() : getProject().getFastaDir());
			}
			set(map);
			if (getLogger().isInfoEnabled()) {
				getLogger().info("Transformed fastq / fasta map: " + map);
			}
		} else {
			set(inputGoal.get());
		}
	}

	protected void fillMap(Map<String, StreamingResourceStream> map, File dir) {
		for (String key : inputGoal.get().keySet()) {
			StreamingResourceListStream list = new StreamingResourceListStream();
			for (StreamingResource resource : inputGoal.get().get(key)) {
				if (resource instanceof StreamingURLResource) {
					StreamingURLResource urlRes = (StreamingURLResource) resource;

					File file = getProject().getOutputFile(dir, null, null, key, null, fastqType ? FileType.FASTQ : FileType.FASTA,
							!urlRes.isNoGZ());
					resource = new StreamingFileResource(file);
				}
				list.getList().add(resource);
			}
			map.put(key, list);
		}
	}
}
