package org.metagene.genestrip;

import java.io.File;
import java.io.IOException;

import org.junit.Test;

public class APITest {
	@Test
	public void testAndDemoAPI() throws IOException {
		File projectDir = getBaseDir();
		
		GSConfig  config = new GSConfig(projectDir);
		File fastq = new File(getBaseDir(), "projects/human_virus/fastq/sample.fastq.gz");
		
		GSProject project = new GSProject(config, "human_virus", 31, null, fastq, null, null);
		GSMaker maker = new GSMaker(project);
		
		maker.getGoal("match").make();
	}
	
	protected File getBaseDir() {
		String projectDir = System.getProperty("project.directory");
		if (projectDir != null) {
			return new File(projectDir, "data");
		}
		String relPath = getClass().getProtectionDomain().getCodeSource().getLocation().getFile();
		return new File(new File(relPath).getParentFile().getParentFile(), "data");
	}	
}
