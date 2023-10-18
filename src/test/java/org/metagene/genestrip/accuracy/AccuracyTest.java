package org.metagene.genestrip.accuracy;

import java.io.File;
import java.io.IOException;

import org.junit.Test;
import org.metagene.genestrip.GSConfig;
import org.metagene.genestrip.GSMaker;
import org.metagene.genestrip.GSProject;

public class AccuracyTest {
	@Test
	public void testAccuracy() throws IOException {
		File testBaseDir =  getTargetDir();
		GSConfig config = new GSConfig(testBaseDir);
		GSProject project = new GSProject(config, "accuracy", 31, null, null, null, null);

		GSMaker maker = new GSMaker(project);
		
		new AccuracyProjectGoal(project, "accproject").make();
		
		AccuracyDataDownloadGoal accDownloadGoal = new AccuracyDataDownloadGoal(project, "accdatadownload", maker.getGoal("setup"));
		FastaTransformGoal fastaTransformGoal = new FastaTransformGoal(project, "fasttransform", accDownloadGoal);

		fastaTransformGoal.make();
	}

	protected File getTargetDir() {
		String projectDir = System.getProperty("project.build.directory");
		if (projectDir != null) {
			return new File(projectDir);
		}
		String relPath = getClass().getProtectionDomain().getCodeSource().getLocation().getFile();
		return new File(relPath).getParentFile();
	}
}
