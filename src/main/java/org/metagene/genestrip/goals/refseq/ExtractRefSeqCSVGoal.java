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
            ps.println("refseq descr; taxid;");
            for (String key : descr2TaxId.keySet()) {
                ps.print(key);
                ps.print(";");
                ps.print(descr2TaxId.get(key));
                ps.println(";");
//                ps.println();
            }
        }
    }
}
