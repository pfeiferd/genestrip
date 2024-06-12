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
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import org.metagene.genestrip.GSGoalKey;
import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.GSProject.FileType;
import org.metagene.genestrip.io.StreamingFileResource;
import org.metagene.genestrip.io.StreamingResource;
import org.metagene.genestrip.io.StreamingURLResource;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;

public class FastqMapTransformGoal extends ObjectGoal<Map<String, List<StreamingResource>>, GSProject> {
	private final ObjectGoal<Map<String, List<StreamingResource>>, GSProject> inputGoal;

	@SafeVarargs
	public FastqMapTransformGoal(GSProject project,
			ObjectGoal<Map<String, List<StreamingResource>>, GSProject> inputGoal, Goal<GSProject>... deps) {
		super(project, GSGoalKey.FASTQ_MAP_TRANSFORM, append(deps, inputGoal));
		this.inputGoal = inputGoal;
	}

	@Override
	public void makeThis() {
		if (getProject().isDownloadFastqs() || getProject().isDownloadFastqsToCommon()) {
			// Linked hash map preserve order of keys as entered.
			Map<String, List<StreamingResource>> map = new LinkedHashMap<String, List<StreamingResource>>();
			if (getProject().isDownloadFastqs()) {
				fillMap(map, getProject().getFastqDir());
			}
			if (getProject().isDownloadFastqsToCommon()) {
				fillMap(map, getProject().getCommon().getFastqDir());
			}
			set(map);
			if (getLogger().isInfoEnabled()) {
				getLogger().info("Transformed fastq map: " + map);
			}
		} else {
			set(inputGoal.get());
		}
	}

	protected void fillMap(Map<String, List<StreamingResource>> map, File dir) {
		for (String key : inputGoal.get().keySet()) {
			List<StreamingResource> list = new ArrayList<StreamingResource>();
			for (StreamingResource resource : inputGoal.get().get(key)) {
				if (resource instanceof StreamingURLResource) {
					StreamingURLResource urlRes = (StreamingURLResource) resource;

					File file = getProject().getOutputFile(dir, null, null, key, null, FileType.FASTQ,
							!urlRes.isNoGZ());
					resource = new StreamingFileResource(file);
				}
				list.add(resource);
			}
			map.put(key, list);
		}
	}
}
