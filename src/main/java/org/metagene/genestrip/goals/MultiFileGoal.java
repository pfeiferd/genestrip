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
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.GSProject.FileType;
import org.metagene.genestrip.io.StreamingResourceStream;
import org.metagene.genestrip.make.FileListGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.GoalKey;
import org.metagene.genestrip.make.ObjectGoal;

public abstract class MultiFileGoal extends FileListGoal<GSProject> {
	protected final ObjectGoal<Map<String, StreamingResourceStream>, GSProject> fastqMapGoal;
	protected final Map<File, StreamingResourceStream> fileToFastqs;
	protected final Map<File, String> fileToKeyMap;

	@SafeVarargs
	public MultiFileGoal(GSProject project, GoalKey key,
			ObjectGoal<Map<String, StreamingResourceStream>, GSProject> fastqMapGoal, Goal<GSProject>... deps) {
		super(project, key, (List<File>) null, append(deps, fastqMapGoal));
		this.fastqMapGoal = fastqMapGoal;
		fileToFastqs = new HashMap<File, StreamingResourceStream>();
		fileToKeyMap = new HashMap<File, String>();
	}

	protected File getSourceDir() {
		return getProject().getFastqDir();
	}

	public Map<File, String> getFileToKeyMap() {
		return fileToKeyMap;
	}

	@Override
	protected void provideFiles() {
		Map<String, StreamingResourceStream> keyToFastqs = fastqMapGoal.get();
		for (String key : keyToFastqs.keySet()) {
			File matchFile = getProject().getOutputFile(getKey().getName(), key, null, getFileType(), isUseGZip());
			addFile(matchFile);
			fileToFastqs.put(matchFile, keyToFastqs.get(key));
			fileToKeyMap.put(matchFile, key);
		}
	}

	protected boolean isUseGZip() {
		return false;
	}

	protected abstract FileType getFileType();
}
