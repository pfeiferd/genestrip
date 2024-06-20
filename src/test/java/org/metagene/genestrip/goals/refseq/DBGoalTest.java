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

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;

import org.junit.BeforeClass;
import org.junit.Test;
import org.metagene.genestrip.APITest;
import org.metagene.genestrip.GSCommon;
import org.metagene.genestrip.GSGoalKey;
import org.metagene.genestrip.GSMaker;
import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.make.FileListGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.make.GoalKey.DefaultGoalKey;
import org.metagene.genestrip.store.Database;

import it.unimi.dsi.fastutil.objects.Object2LongMap;

public class DBGoalTest {
	@BeforeClass
	public static void clearDB() throws IOException {
		GSMaker maker = new GSMaker(createProject());

		maker.getGoal(GSGoalKey.CLEAR).make();
		maker.dumpAll();
	}

	protected static GSProject createProject() throws IOException {
		File baseDir = getBaseDir();

		// Load the system configuration.
		GSCommon config = new GSCommon(baseDir) {
			@Override
			public File getCommonDir() {
				return new File(APITest.getBaseDir(), "common");
			}
		};

		return new GSProject(config, "dengue1");
	}

	protected static File getBaseDir() {
		return new File(getTargetDir(), "data");
	}

	protected static File getTargetDir() {
		String projectDir = System.getProperty("project.build.directory");
		if (projectDir != null) {
			return new File(projectDir);
		}
		String relPath = DBGoalTest.class.getProtectionDomain().getCodeSource().getLocation().getFile();
		return new File(relPath).getParentFile();
	}
	
	@Test
	public void testUpdate() throws IOException {
		GSProject project = createProject();

		new Dengue1ProjectGoal(project).make();
		
		GSMaker maker = new GSMaker(project);

		@SuppressWarnings("unchecked")
		ObjectGoal<Database, GSProject> fillGoal = (ObjectGoal<Database, GSProject>) maker.getGoal(GSGoalKey.FILL_DB);
		Object2LongMap<String> stats = fillGoal.get().getStats();
		long totalKmers = stats.getLong(null);
		assertTrue(totalKmers > 0);

		long kmers = stats.getLong("11053");
		assertEquals(totalKmers, kmers);

		@SuppressWarnings("unchecked")
		ObjectGoal<Database, GSProject> dbGoal = (ObjectGoal<Database, GSProject>) maker.getGoal(GSGoalKey.UPDATE_DB);
		Object2LongMap<String> stats2 = dbGoal.get().getStats();
		
		long kmers2 = stats2.getLong("11053");
		assertEquals(0, kmers2);
		long totalKmers2 = stats2.getLong(null);
		assertEquals(totalKmers, totalKmers2);
		long kmers3 = stats2.getLong("1");
		assertEquals(totalKmers, kmers3);
		
		// Clean up memory and threads.
		maker.dumpAll();
	}
	
	private static class Dengue1ProjectGoal extends FileListGoal<GSProject> {
		@SafeVarargs
		public Dengue1ProjectGoal(GSProject project, Goal<GSProject>... dependencies) {
			super(project, new DefaultGoalKey("dengue1"), (List<File>) null, dependencies);
			addFile(new File(project.getProjectDir(), "taxids.txt"));
			addFile(new File(project.getProjectDir(), "categories.txt"));
			addFile(new File(project.getProjectDir(), "additional.txt"));
			addFile(new File(project.getProjectDir(), "fasta/dengue1.fasta"));
		}

		@Override
		protected void makeFile(File file) throws IOException {
			if (!getProject().getCommon().getBaseDir().exists()) {
				getProject().getCommon().getBaseDir().mkdir();
			}
			if (!getProject().getProjectsDir().exists()) {
				getProject().getProjectsDir().mkdir();
			}
			if (!getProject().getProjectDir().exists()) {
				getProject().getProjectDir().mkdir();
			}
			if (!getProject().getFastaDir().exists()) {
				getProject().getFastaDir().mkdir();
			}
			if (!file.exists()) {
				ClassLoader classLoader = getClass().getClassLoader();
				File orgFile = new File(classLoader.getResource("projects/dengue1/" + file.getName()).getFile());
				Files.copy(orgFile.toPath(), file.toPath());
			}
		}
	}	
}
