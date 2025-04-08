package org.metagene.genestrip.goals;

import org.metagene.genestrip.GSConfigKey;
import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.fasta.AbstractFastaReader;
import org.metagene.genestrip.io.StreamProvider;
import org.metagene.genestrip.io.StreamingResource;
import org.metagene.genestrip.io.StreamingResourceStream;
import org.metagene.genestrip.make.FileListGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.GoalKey;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.util.ByteArrayUtil;

import java.io.*;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class Fasta2FastqGoal extends FileListGoal<GSProject> {
    private final ObjectGoal<Map<String, StreamingResourceStream>, GSProject> fastaMapGoal;
    private final Map<File, String> fileToKeyMap;

    public Fasta2FastqGoal(GSProject project, GoalKey key, ObjectGoal<Map<String, StreamingResourceStream>, GSProject> fastaMapGoal, Goal<GSProject>... deps) {
        super(project, key, (List<File>) null, append(deps, fastaMapGoal));

        this.fastaMapGoal = fastaMapGoal;
        this.fileToKeyMap = new HashMap<>();
    }

    @Override
    protected void provideFiles() {
        for (String key : fastaMapGoal.get().keySet()) {
            File fastqFile = getProject().getOutputFile(getKey().getName(), key, null, GSProject.FileType.FASTQ, true);
            addFile(fastqFile);
            fileToKeyMap.put(fastqFile, key);
        }
    }


    @Override
    protected void makeFile(File file) throws IOException {
        try {
            StreamingResourceStream fastas = fastaMapGoal.get().get(fileToKeyMap.get(file));
            try (PrintStream out = new PrintStream(StreamProvider.getOutputStreamForFile(file))) {
                FastqWriter writer = new FastqWriter(out, 65535);
                for (StreamingResource rs : fastas) {
                    writer.readFasta(rs.openStream().getInputStream());
                }
            }
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    protected static class FastqWriter extends AbstractFastaReader {
        private int dataSize;
        private final PrintStream out;

        public FastqWriter(PrintStream out, int bufferSize) {
            super(bufferSize);
            this.out = out;
        }

        @Override
        protected void infoLine() {
            out.print("@");
            ByteArrayUtil.println(target, 0, size - 1, out);
        }

        @Override
        protected void startRegion() {
            dataSize = 0;
        }

        @Override
        protected void dataLine() {
            ByteArrayUtil.print(target, 0, size - 1, out);
            dataSize += size - 1;
        }

        @Override
        protected void endRegion() {
            out.println();
            out.println("+");
            for (int i = 0; i < dataSize; i++) {
                out.print('~');
            }
            out.println();
        }
    }
}
