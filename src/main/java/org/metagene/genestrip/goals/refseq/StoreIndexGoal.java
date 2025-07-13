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
package org.metagene.genestrip.goals.refseq;

import java.io.File;
import java.io.IOException;

import org.metagene.genestrip.GSGoalKey;
import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.GSProject.FileType;
import org.metagene.genestrip.bloom.AbstractKMerBloomFilter;
import org.metagene.genestrip.make.FileListGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;

public class StoreIndexGoal extends FileListGoal<GSProject> {
	private final ObjectGoal<AbstractKMerBloomFilter, GSProject> indexGoal;

	@SafeVarargs
	public StoreIndexGoal(GSProject project, ObjectGoal<AbstractKMerBloomFilter, GSProject> indexGoal,
			Goal<GSProject>... deps) {
		super(project, GSGoalKey.INDEX, project.getOutputFile(GSGoalKey.INDEX.getName(), FileType.FILTER, true),
				Goal.append(deps, indexGoal));
		this.indexGoal = indexGoal;
	}

	@Override
	protected void makeFile(File indexFile) {
		try {
			AbstractKMerBloomFilter index = indexGoal.get();
			if (getLogger().isInfoEnabled()) {
				getLogger().info("Saving index " + indexFile + " ...");
			}
			index.save(indexFile);
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}
}