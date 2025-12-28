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

import me.tongfei.progressbar.ProgressBar;
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
import org.metagene.genestrip.util.progressbar.GSProgressBarCreator;

import java.io.*;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class Fasta2FastqGoal<P extends GSProject> extends FileListGoal<P> {
    private final ObjectGoal<Map<String, StreamingResourceStream>, P> fastaMapGoal;
    private final Map<File, String> fileToKeyMap;

    public Fasta2FastqGoal(P project, GoalKey key, ObjectGoal<Map<String, StreamingResourceStream>, P> fastaMapGoal, Goal<P>... deps) {
        super(project, key, (List<File>) null, append(deps, fastaMapGoal));

        this.fastaMapGoal = fastaMapGoal;
        this.fileToKeyMap = new HashMap<>();
    }

    @Override
    protected void provideFiles() {
        for (String key : fastaMapGoal.get().keySet()) {
            File fastqFile = getProject().getOutputFile(getKey().getName(), key, null, GSProject.GSFileType.FASTQ, true);
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
                    try (StreamingResource.StreamAccess byteCountAccess = rs.openStream()) {
                        try (ProgressBar pb = isProgressBar() ? GSProgressBarCreator.newGSProgressBar(getKey().getName(), byteCountAccess, null, true) : null) {
                            writer.readFasta(byteCountAccess.getInputStream());
                        }
                    }
                }
            }
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    protected boolean isProgressBar() {
        return getProject().booleanConfigValue(GSConfigKey.PROGRESS_BAR);
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
