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

import org.apache.commons.io.FileUtils;
import org.junit.Test;
import org.metagene.genestrip.*;
import org.metagene.genestrip.fastq.AbstractFastqReader;
import org.metagene.genestrip.fastq.AbstractLoggingFastqStreamer;
import org.metagene.genestrip.goals.refseq.DBGoalTest;
import org.metagene.genestrip.io.StreamProvider;
import org.metagene.genestrip.io.StreamingResourceStream;

import java.io.*;

import static org.junit.Assert.assertTrue;

public class Fasta2FastqGoalTest extends DBGoalTest {
    @Override
    public void testUpdate() throws IOException {
        // Just to avoid running the test from the superclass ...
    }

    @Override
    public void testKrakenOutput() {
        // Just to avoid running the test from the superclass ...
    }

    @Test
    public void testFasta2Fast2Goal() throws IOException {
        GSCommon config = new GSCommon(getBaseDir());
        GSProject project = new GSProject(config, "dengue1", null, new String[] { "fasta2fastqtest.fasta" }, null, null, null, null,
                null, null, null, false);
        project.initConfigParam(GSConfigKey.GZIP_FASTQ_OUTPUT, false);

        int k = project.intConfigValue(GSConfigKey.KMER_SIZE);

        createProjectGoal(project).make();

        GSMaker maker = new GSMaker(project);
        Fasta2FastqGoal goal = (Fasta2FastqGoal) maker.getGoal(GSGoalKey.FASTA2FASTQ);
        goal.cleanThis();
        goal.make();

        File testOut = new File(getBaseDir(), "test_out.fastq");

        try (OutputStream ps = new FileOutputStream(testOut)) {
            StreamingResourceStream streamingResourceStream = goal.getFastasForFile(goal.getFile());
            new AbstractLoggingFastqStreamer(k, project.intConfigValue(GSConfigKey.INITIAL_READ_SIZE_BYTES), 0, new DefaultExecutionContext(null, 0, 1),
                    true) {
                @Override
                protected void nextEntry(AbstractFastqReader.ReadEntry readStruct, int threadIndex) throws IOException {
                    readStruct.write(ps);
                }
            }.processFastqStreams(streamingResourceStream);
        }

        assertTrue(FileUtils.contentEquals(testOut, goal.getFile()));

        // Read the output once via fastq reader as a check for structural consistency.
        new AbstractFastqReader(k, project.intConfigValue(GSConfigKey.INITIAL_READ_SIZE_BYTES), 0, new DefaultExecutionContext(null, 0, 1),
                true) {
            @Override
            protected void nextEntry(ReadEntry readStruct, int threadIndex) throws IOException {
            }

            @Override
            public void readFastq(InputStream inputStream, boolean fasta) throws IOException {
                super.readFastq(inputStream, fasta);
            }
        }.readFastq(StreamProvider.getInputStreamForFile(goal.getFile()), false);

        maker.dumpAll();
    }
}
