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

import me.tongfei.progressbar.DelegatingProgressBarConsumer;
import me.tongfei.progressbar.ProgressBar;
import me.tongfei.progressbar.ProgressBarBuilder;
import me.tongfei.progressbar.ProgressBarStyle;
import org.metagene.genestrip.GSConfigKey;
import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.GoalKey;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.refseq.AbstractRefSeqFastaReader;
import org.metagene.genestrip.refseq.RefSeqCategory;
import org.metagene.genestrip.tax.TaxTree;
import org.metagene.genestrip.util.GSLogFactory;
import org.metagene.genestrip.util.progressbar.GSProgressBarCreator;

import java.io.File;
import java.io.IOException;
import java.time.temporal.ChronoUnit;
import java.util.List;
import java.util.Map;
import java.util.Set;

public abstract class FastaReaderGoal<T> extends ObjectGoal<T, GSProject> {
    protected final ObjectGoal<Set<RefSeqCategory>, GSProject> categoriesGoal;
    protected final ObjectGoal<Set<TaxTree.TaxIdNode>, GSProject> taxNodesGoal;
    protected final RefSeqFnaFilesDownloadGoal fnaFilesGoal;
    protected final ObjectGoal<Map<File, TaxTree.TaxIdNode>, GSProject> additionalGoal;

    public FastaReaderGoal(GSProject project, GoalKey key, ObjectGoal<Set<RefSeqCategory>, GSProject> categoriesGoal,
                           ObjectGoal<Set<TaxTree.TaxIdNode>, GSProject> taxNodesGoal, RefSeqFnaFilesDownloadGoal fnaFilesGoal,
                           ObjectGoal<Map<File, TaxTree.TaxIdNode>, GSProject> additionalGoal, Goal<GSProject>... dependencies) {
        super(project, key, Goal.append(dependencies, taxNodesGoal, fnaFilesGoal, additionalGoal));
        this.categoriesGoal = categoriesGoal;
        this.taxNodesGoal = taxNodesGoal;
        this.fnaFilesGoal = fnaFilesGoal;
        this.additionalGoal = additionalGoal;
    }

    protected void readFastas(AbstractRefSeqFastaReader fastaReader) throws IOException {
        boolean refSeqDB = booleanConfigValue(GSConfigKey.REF_SEQ_DB);
        int sumFiles = 0;
        List<File> refSeqFiles = null;
        if (refSeqDB) {
            fnaFilesGoal.make();
            refSeqFiles = fnaFilesGoal.getFiles();
            sumFiles += refSeqFiles.size();
        }
        Map<File, TaxTree.TaxIdNode> additionalMap = additionalGoal.get();
        sumFiles += additionalMap.size();
        try (ProgressBar pb = createProgressBar(sumFiles)) {
            if (refSeqFiles != null) {
                for (File fnaFile : refSeqFiles) {
                    RefSeqCategory cat = fnaFilesGoal.getCategoryForFile(fnaFile);
                    if (categoriesGoal.get().contains(cat)) {
                        fastaReader.readFasta(fnaFile);
                        if (pb != null) {
                            pb.step();
                        }
                    }
                }
            }
            for (File additionalFasta : additionalMap.keySet()) {
                TaxTree.TaxIdNode node = additionalMap.get(additionalFasta);
                if (taxNodesGoal.get().contains(node)) {
                    fastaReader.ignoreAccessionMap(node);
                    fastaReader.readFasta(additionalFasta);
                    if (pb != null) {
                        pb.step();
                    }
                }
            }
        }
    }

    protected ProgressBar createProgressBar(int max) {
        return booleanConfigValue(GSConfigKey.PROGRESS_BAR) ?
                GSProgressBarCreator.newGSProgressBar(getKey().getName(), 60000, " files", null, getLogger()) :
                null;
    }

    @Override
    public boolean isWeakDependency(Goal<GSProject> toGoal) {
        if (toGoal == fnaFilesGoal) {
            return true;
        }
        return super.isWeakDependency(toGoal);
    }
}
