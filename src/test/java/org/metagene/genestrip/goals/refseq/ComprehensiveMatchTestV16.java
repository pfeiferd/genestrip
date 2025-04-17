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
import org.junit.Ignore;
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

@Ignore
public class ComprehensiveMatchTestV16 extends ComprehensiveMatchTest {
    @BeforeClass()
    public static void clearDB() throws IOException {
        GSProject project = createProject("viral16", null);
        GSMaker maker = new GSMaker(project);
        maker.getGoal(GSGoalKey.CLEAR).make();
        maker.dumpAll();
    }

    protected String getProjectName() {
        return "viral16";
    }

    @Override
    protected FileListGoal<GSProject> createProjectGoal(GSProject project) {
        return new ViralProjectGoal16(project);
    }

    @Override
    protected GSProject createTestProject(String csvFile1) throws IOException {
        GSProject project = super.createTestProject(csvFile1);
        project.initConfigParam(GSConfigKey.KMER_SIZE, 16);
        return project;
    }

    public class ViralProjectGoal16 extends ViralProjectGoal {
        @SafeVarargs
        public ViralProjectGoal16(GSProject project, Goal<GSProject>... dependencies) {
            super(project, dependencies);
        }

        @Override
        protected String getResFolderStr(String fileName) {
            return "test.fasta.gz".equals(fileName) ? "viral" : getProjectName();
        }
    }
}
