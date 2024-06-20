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
package org.metagene.genestrip;

import java.io.File;
import java.io.IOException;

import org.junit.BeforeClass;
import org.junit.Test;
import org.metagene.genestrip.goals.GSFileDownloadGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.GoalKey;

public class GoalsTest {
	@BeforeClass
	public static void clearDB() throws IOException {
		APITest.clearDB();
	}

	@Test
	public void testGoals() throws IOException {
		File baseDir = APITest.getBaseDir();
		GSCommon config = new GSCommon(baseDir);
		File fastq = new File(APITest.getBaseDir(), "projects/human_virus/fastq/sample.fastq.gz");
		GSProject project = new GSProject(config, "human_virus", null, new String[] { fastq.getCanonicalPath() });
		GSMaker maker = new GSMaker(project);
		for (Goal<GSProject> goal : maker.getGoals()) {
			GoalKey key = goal.getKey();
			if (!GSGoalKey.KRAKENRES.equals(key)) {
				// Don't download NCBI stuff over and over again via tests:
				if (!GSGoalKey.COMMON_SETUP.equals(key) && !(goal instanceof GSFileDownloadGoal)) {
					goal.cleanThis();
				}
				goal.make();
			}
		}
		maker.dumpAll();
	}
}
