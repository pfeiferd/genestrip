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
import org.metagene.genestrip.GSProject.GSFileType;
import org.metagene.genestrip.io.StreamingResourceStream;
import org.metagene.genestrip.make.FileListGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.GoalKey;
import org.metagene.genestrip.make.ObjectGoal;

/**
 * Base class for goals that produce one output file per fastq/fasta map key, keeping the association from each
 * output file back to its input resources.
 *
 * @param <P> the project type
 */
public abstract class MultiFileGoal<P extends GSProject> extends FileListGoal<P> {
	/** The goal supplying the per-key fastq/fasta input streams. */
	protected final ObjectGoal<Map<String, StreamingResourceStream>, P> fastqMapGoal;
	/** Maps each produced output file back to its input resources. */
	protected final Map<File, StreamingResourceStream> fileToFastqs;

	/**
	 * Creates the goal wiring the fastq map goal whose keys drive the produced output files.
	 *
	 * @param project the project this goal belongs to
	 * @param key the goal key identifying this goal
	 * @param fastqMapGoal the goal supplying the per-key fastq/fasta input streams
	 * @param deps additional goals this goal depends on
	 */
	@SafeVarargs
	public MultiFileGoal(P project, GoalKey key,
			ObjectGoal<Map<String, StreamingResourceStream>, P> fastqMapGoal, Goal<P>... deps) {
		super(project, key, (List<File>) null, append(deps, fastqMapGoal));
		this.fastqMapGoal = fastqMapGoal;
		fileToFastqs = new HashMap<>();
	}

	/*
	public Map<File, String> getFileToKeyMap() {
		return fileToKeyMap;
	}
	*/

	@Override
	protected void provideFiles() {
		Map<String, StreamingResourceStream> keyToFastqs = fastqMapGoal.get();
		for (String key : keyToFastqs.keySet()) {
			File matchFile = getProject().getOutputFile(getKey().getName(), key, null, getFileType(), isUseGZip());
			addFile(matchFile);
			fileToFastqs.put(matchFile, keyToFastqs.get(key));
			enterFileAndKey(matchFile, key);
		}
	}

	/**
	 * Hook called for each provided output file together with its map key; does nothing by default.
	 *
	 * @param file the produced output file
	 * @param key the map key the file corresponds to
	 */
	protected void enterFileAndKey(File file, String key) {}

	/**
	 * Returns whether produced output files should be GZIP-compressed.
	 *
	 * @return whether output files should be GZIP-compressed
	 */
	protected boolean isUseGZip() {
		return false;
	}

	/**
	 * Returns the file type of the produced output files.
	 *
	 * @return the output file type
	 */
	protected abstract GSFileType getFileType();
}
