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
package org.metagene.genestrip.finertree.goals;

import org.metagene.genestrip.finertree.FTConfigKey;
import org.metagene.genestrip.finertree.FTGoalKey;
import org.metagene.genestrip.finertree.FTProject;
import org.metagene.genestrip.make.FileGoal;
import org.metagene.genestrip.make.Goal;

import java.io.*;
import java.nio.file.Files;
import java.util.*;

public class AllInOneLaTeXGoal<P extends FTProject>  extends FileGoal<P> {
    private final DengrogramLaTeXGoal<P> dengrogramLaTeXGoal;
    private final Map<File, List<File>> outToInFiles;

    public AllInOneLaTeXGoal(P project, DengrogramLaTeXGoal<P> dengrogramLaTeXGoal, Goal<P>... deps) {
        super(project, FTGoalKey.ALLINONE_LATEX, append(deps, dengrogramLaTeXGoal));
        this.dengrogramLaTeXGoal = dengrogramLaTeXGoal;
        this.outToInFiles = new HashMap<>();
    }

    @Override
    public boolean isCleaned() {
        // Actual check would require make of other goals - not doable.
        return true;
    }

    @Override
    public List<File> getFiles() {
        List<File> files = dengrogramLaTeXGoal.getFiles();
        int total = files.size();
        int stepSize = intConfigValue(FTConfigKey.ALLINONE_CHUNK_SIZE);
        List<File> res = new ArrayList<>();
        for (int i = 0; i < total; i += stepSize) {
            File outputFile = getProject().getOutputFile(getKey().getName(), Integer.toString(i), null, FTProject.FTFileType.TEX, false);
            res.add(outputFile);
            List<File> inFiles = new ArrayList<>();
            for (int j = i; j < i + stepSize && j < total; j++) {
                inFiles.add(files.get(j));
            }
            outToInFiles.put(outputFile, inFiles);
        }
        return res;
    }

    @Override
    protected void makeFile(File file) throws IOException {
        try (FileOutputStream out = new FileOutputStream(file)) {
            try (PrintStream pout = new PrintStream(out)) {
                //pout.println("\\documentclass[border=0]{standalone}");
                pout.println("\\documentclass{article}");
                pout.println("\\usepackage[paperwidth=21cm,paperheight=200cm,margin=1cm]{geometry}");
                pout.println("\\usepackage{tikz}");
                pout.println("\\begin{document}");
                //pout.println("\\begin{minipage}{21cm}");
                pout.flush();
                for (File latexFile : outToInFiles.get(file)) {
                    Files.copy(latexFile.toPath(), out);
                }
                out.flush();
                //pout.println("\\end{minipage}");
                pout.println("\\end{document}");
            }
        }
    }
}
