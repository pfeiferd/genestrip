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

import org.apache.commons.io.FileUtils;
import org.junit.BeforeClass;
import org.junit.Test;
import org.metagene.genestrip.*;
import org.metagene.genestrip.make.FileListGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.GoalKey;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;

import static org.junit.Assert.assertTrue;

public class ComprehensiveMatchTest extends DBGoalTest {
    @BeforeClass()
    public static void clearDB() throws IOException {
        GSProject project = createProject("viral", null);
        GSMaker maker = new GSMaker(project);
        maker.getGoal(GSGoalKey.CLEAR).make();
        maker.dumpAll();
    }

    protected GSProject createTestProject(String csvFile1) throws IOException {
        return createProject(getProjectName(), csvFile1);
    }

    protected String getProjectName() {
        return "viral";
    }

    @Override
    protected FileListGoal<GSProject> createProjectGoal(GSProject project) {
        return new ViralProjectGoal(project);
    }

    @Override
    public void testUpdate() throws IOException {
        // Just to avoid running the test from the superclass ...
    }

    @Test
    public void testKrakenOutput() throws IOException {
        GSProject project = createTestProject("viral_test_fasta.txt");
        createProjectGoal(project).make();

        GSMaker maker = new GSMaker(project);
        maker.getGoal(GSGoalKey.FASTA2FASTQ).make();
        maker.dumpAll();

        project = createTestProject(null);
        maker = new GSMaker(project);
        maker.getGoal(GSGoalKey.FASTA2FASTQ).make();
        maker.getGoal(GSGoalKey.DB).make();
        maker.dumpAll();

        for (int i = 0; i < 2; i++) {
            project = createTestProject("viral_test_fastq.txt");
            project.initConfigParam(GSConfigKey.THREADS, 0);
            project.initConfigParam(GSConfigKey.WRITE_KRAKEN_STYLE_OUT, true);
            project.initConfigParam(GSConfigKey.USE_INLINED, i == 0);
            maker = new GSMaker(project);
            maker.match(false, "test", new File(project.getFastqDir(), getProjectName() +  "_fasta2fastq_test.fastq.gz").toString());
            File file1 = new File(project.getKrakenOutDir(), "test_ks.out");
            File file2 = new File(project.getKrakenOutDir(), getProjectName() + "_matchres_test.out");
            System.out.println("file1: " + file1);
            System.out.println("file2: " + file2);
            assertTrue(FileUtils.contentEquals(file1, file2));
            // Clean up memory and threads.
            maker.dumpAll();
        }
    }

    public class ViralProjectGoal extends FileListGoal<GSProject> {
        @SafeVarargs
        public ViralProjectGoal(GSProject project, Goal<GSProject>... dependencies) {
            super(project, new GoalKey.DefaultGoalKey(getProjectName()), (List<File>) null, dependencies);
            addFile(new File(project.getProjectDir(), "taxids.txt"));
            addFile(new File(project.getProjectDir(), "categories.txt"));
            addFile(new File(project.getProjectDir(), "fasta/test.fasta.gz"));
            addFile(new File(project.getProjectDir(), "fastq/viral_test_fastq.txt"));
            addFile(new File(project.getProjectDir(), "fasta/viral_test_fasta.txt"));
            addFile(new File(project.getProjectDir(), "krakenout/test_ks.out"));
            addFile(new File(project.getProjectDir(), "krakenout/test_gs.out"));
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
            if (!getProject().getFastqDir().exists()) {
                getProject().getFastqDir().mkdir();
            }
            if (!getProject().getKrakenOutDir().exists()) {
                getProject().getKrakenOutDir().mkdir();
            }
            if (!file.exists()) {
                ClassLoader classLoader = getClass().getClassLoader();
                File orgFile = new File(classLoader.getResource("projects/" + getResFolderStr(file.getName()) + "/" + file.getName()).getFile());
                Files.copy(orgFile.toPath(), file.toPath());
            }
        }

        protected String getResFolderStr(String fileName) {
            return getProjectName();
        }
    }
}
