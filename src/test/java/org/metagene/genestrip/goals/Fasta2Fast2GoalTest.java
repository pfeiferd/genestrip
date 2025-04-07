package org.metagene.genestrip.goals;

import org.junit.Test;
import org.metagene.genestrip.*;
import org.metagene.genestrip.fastq.AbstractFastqReader;
import org.metagene.genestrip.goals.refseq.DBGoalTest;
import org.metagene.genestrip.io.StreamProvider;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;

public class Fasta2Fast2GoalTest extends DBGoalTest {
    @Test
    public void testFasta2Fast2Goal() throws IOException {
        File baseDir = APITest.getBaseDir();

        GSCommon config = new GSCommon(DBGoalTest.getBaseDir());
        GSProject project = new GSProject(config, "dengue1", null, new String[] { "fasta2fastqtest.fasta" }, null, null, null, false, null,
                null, null, null, false);

        GSMaker maker = new GSMaker(project);

        Fasta2FastqGoal goal = (Fasta2FastqGoal) maker.getGoal(GSGoalKey.FASTA2FASTQ);
        goal.cleanThis();
        new Dengue1ProjectGoal(project).make();
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
    }
}
