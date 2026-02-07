package org.metagene.genestrip.goals.refseq;

import org.metagene.genestrip.ExecutionContext;
import org.metagene.genestrip.GSConfigKey;
import org.metagene.genestrip.GSGoalKey;
import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.io.StreamProvider;
import org.metagene.genestrip.make.FileGoal;
import org.metagene.genestrip.make.FileListGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.refseq.AbstractRefSeqFastaReader;
import org.metagene.genestrip.refseq.AccessionMap;
import org.metagene.genestrip.refseq.RefSeqCategory;
import org.metagene.genestrip.tax.Rank;
import org.metagene.genestrip.tax.TaxTree;
import org.metagene.genestrip.util.ByteArrayUtil;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.*;

public class ExtractRefSeqCSVGoal<P extends GSProject> extends FileListGoal<P> {
    private final ObjectGoal<Map<String, String>, P> extractFastaGoal;

    @SafeVarargs
    public ExtractRefSeqCSVGoal(P project, ObjectGoal<Map<String, String>, P> extractFastaGoal, Goal<P>... deps) {
        super(project, GSGoalKey.EXTRACT_REFSEQ_CSV, project.getOutputFile(GSGoalKey.EXTRACT_REFSEQ_CSV.getName(), GSProject.GSFileType.CSV, false), Goal.append(deps, extractFastaGoal));
        this.extractFastaGoal = extractFastaGoal;
    }

    @Override
    protected void makeFile(File file) throws IOException {
        Map<String, String> descr2TaxId = extractFastaGoal.get();
        try (PrintStream ps = new PrintStream(StreamProvider.getOutputStreamForFile(file))) {
            //ps.println("refseq descr; taxid;");
            for (String key : descr2TaxId.keySet()) {
                ps.print(key);
                ps.print("\t");
//                    ps.print(";");
                ps.print(descr2TaxId.get(key));
//                    ps.println(";");
                ps.println();
            }
        }
    }
}
