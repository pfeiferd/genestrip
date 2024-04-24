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
