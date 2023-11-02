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
		addFile(new File(project.getProjectDir(), "project.properties"));
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
