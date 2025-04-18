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
package org.metagene.genestrip.goals;

import org.junit.Ignore;
import org.junit.Test;
import org.metagene.genestrip.*;
import org.metagene.genestrip.fastq.AbstractFastqReader;
import org.metagene.genestrip.goals.refseq.DBGoalTest;
import org.metagene.genestrip.io.StreamProvider;

import java.io.IOException;
import java.io.InputStream;

public class Fasta2Fast2GoalTest extends DBGoalTest {
    @Override
    public void testUpdate() throws IOException {
        // Just to avoid running the test from the superclass ...
    }

    @Override
    public void testKrakenOutput() throws IOException {
        // Just to avoid running the test from the superclass ...
    }

    @Test
    public void testFasta2Fast2Goal() throws IOException {
        GSCommon config = new GSCommon(getBaseDir());
        GSProject project = new GSProject(config, "dengue1", null, new String[] { "fasta2fastqtest.fasta" }, null, null, null, null,
                null, null, null, false);

        createProjectGoal(project).make();

        GSMaker maker = new GSMaker(project);
        Fasta2FastqGoal goal = (Fasta2FastqGoal) maker.getGoal(GSGoalKey.FASTA2FASTQ);
        goal.cleanThis();
        goal.make();

        new AbstractFastqReader(31, project.intConfigValue(GSConfigKey.INITIAL_READ_SIZE_BYTES), 0, new DefaultExecutionContext(0, 1),
                true) {
            @Override
            protected void nextEntry(ReadEntry readStruct, int threadIndex) throws IOException {
            }

            @Override
            public void readFastq(InputStream inputStream) throws IOException {
                super.readFastq(inputStream);
            }
        }.readFastq(StreamProvider.getInputStreamForFile(goal.getFile()));

        maker.dumpAll();
    }
}
