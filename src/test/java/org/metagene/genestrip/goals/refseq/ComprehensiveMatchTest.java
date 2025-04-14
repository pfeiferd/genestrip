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
    @Override
    protected FileListGoal<GSProject> createProjectGoal(GSProject project) {
        return new ViralProjectGoal(project);
    }

    @Override
    public void testUpdate() throws IOException {
        // Just to avoid running the test from the superclass...
    }

    protected static GSProject createProject(String csvFile1) throws IOException {
        File baseDir = getBaseDir();

        // Load the system configuration.
        GSCommon config = new GSCommon(baseDir) {
            @Override
            public File getCommonDir() {
                return new File(APITest.getBaseDir(), "common");
            }
        };

        return new GSProject(config, "viral", null, null, csvFile1, null, null, false, null,
                null, null, null, false);
    }


    @Test
    public void testKrakenOutput() throws IOException {
        GSProject project = createProject("viral_test_fasta.txt");
        createProjectGoal(project).make();

        GSMaker maker = new GSMaker(project);
        maker.getGoal(GSGoalKey.FASTA2FASTQ).make();
        maker.dumpAll();

        project = createProject("viral_test_fastq.txt");
        project.initConfigParam(GSConfigKey.WRITED_KRAKEN_STYLE_OUT, true);

        createProjectGoal(project).make();

        maker = new GSMaker(project);
        maker.getGoal(GSGoalKey.FASTA2FASTQ).make();

        maker.match(false, "test", new File(project.getFastqDir(), "viral_fasta2fastq_test.fastq.gz").toString());
        File file1 = new File(project.getKrakenOutDir(), "test.out");
        File file2 = new File(project.getKrakenOutDir(), "viral_test.out");
        assertTrue(FileUtils.contentEquals(file1, file2));
        // Clean up memory and threads.
        maker.dumpAll();
    }


    public static class ViralProjectGoal extends FileListGoal<GSProject> {
        @SafeVarargs
        public ViralProjectGoal(GSProject project, Goal<GSProject>... dependencies) {
            super(project, new GoalKey.DefaultGoalKey("viral"), (List<File>) null, dependencies);
            addFile(new File(project.getProjectDir(), "taxids.txt"));
            addFile(new File(project.getProjectDir(), "categories.txt"));
            addFile(new File(project.getProjectDir(), "fasta/test.fasta.gz"));
            addFile(new File(project.getProjectDir(), "fastq/viral_test_fastq.txt"));
            addFile(new File(project.getProjectDir(), "fasta/viral_test_fasta.txt"));
            addFile(new File(project.getProjectDir(), "krakenout/test.out"));
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
                File orgFile = new File(classLoader.getResource("projects/viral/" + file.getName()).getFile());
                Files.copy(orgFile.toPath(), file.toPath());
            }
        }
    }
}
