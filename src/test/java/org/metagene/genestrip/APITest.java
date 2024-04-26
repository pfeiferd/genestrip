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

import org.junit.Test;
import org.metagene.genestrip.make.Goal;

public class APITest {
	@Test
	public void testAndDemoAPI() throws IOException {
		// Get the base directory for Genestrip. (Depends on your setup.)
		File baseDir = getBaseDir();

		// Load the system configuration.
		GSConfig config = new GSConfig(baseDir);

		// Get the fastq file to run a match on. (In this case, a small sample file
		// provied as part of the release.)
		File fastq = new File(getBaseDir(), "projects/human_virus/fastq/sample.fastq.gz");

		// Create the 'human_virus' project (whose config files are part of the
		// release).
		GSProject project = new GSProject(config, "human_virus", 31, null, fastq, null, null, null);
		GSMaker maker = new GSMaker(project);
		// Run the 'match' goal. This may trigger other goals such as the 'db' goal to
		// first create the 'human_virus' database.
		Goal<GSProject> goal = maker.getGoal(GSMaker.UserGoal.MATCH);
		goal.cleanThis();
		goal.make();
	}

	// To find the base directory for any project file provided as part of the
	// release.
	// (This may have to be done differently for you own projects.)
	public static File getBaseDir() {
		String projectDir = System.getProperty("project.directory");
		if (projectDir != null) {
			return new File(projectDir, "data");
		}
		String relPath = APITest.class.getProtectionDomain().getCodeSource().getLocation().getFile();
		return new File(new File(relPath).getParentFile().getParentFile(), "data");
	}
}
