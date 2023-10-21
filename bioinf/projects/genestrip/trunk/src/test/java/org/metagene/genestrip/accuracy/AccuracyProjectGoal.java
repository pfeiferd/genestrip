package org.metagene.genestrip.accuracy;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;

import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.make.FileListGoal;
import org.metagene.genestrip.make.Goal;

public class AccuracyProjectGoal extends FileListGoal<GSProject> {
	@SafeVarargs
	public AccuracyProjectGoal(GSProject project, String name, Goal<GSProject>... dependencies) {
		super(project, name, (List<File>) null, dependencies);
// We generate the taxids.txt file instead because we need the taxid on the genus level.		
//		addFile(new File(project.getProjectDir(), "taxids.txt"));
		addFile(new File(project.getProjectDir(), "categories.txt"));
		addFile(new File(project.getProjectDir(), "multimatch.txt"));
	}

	@Override
	protected void makeFile(File file) throws IOException {
		if (!getProject().getProjectsDir().exists()) {
			getProject().getProjectsDir().mkdir();
		}
		if (!getProject().getProjectDir().exists()) {
			getProject().getProjectDir().mkdir();
		}
		if (!file.exists()) {
			ClassLoader classLoader = getClass().getClassLoader();
			File orgFile = new File(classLoader.getResource("projects/accuracy/" + file.getName()).getFile());
			Files.copy(orgFile.toPath(), file.toPath());
		}
	}
}
