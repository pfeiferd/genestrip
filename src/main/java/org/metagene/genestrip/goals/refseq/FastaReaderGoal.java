package org.metagene.genestrip.goals.refseq;

import org.metagene.genestrip.GSConfigKey;
import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.GoalKey;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.refseq.AbstractRefSeqFastaReader;
import org.metagene.genestrip.refseq.RefSeqCategory;
import org.metagene.genestrip.tax.TaxTree;

import java.io.File;
import java.io.IOException;
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
        if (refSeqDB) {
            fnaFilesGoal.make();
            for (File fnaFile : fnaFilesGoal.getFiles()) {
                RefSeqCategory cat = fnaFilesGoal.getCategoryForFile(fnaFile);
                if (categoriesGoal.get().contains(cat)) {
                    fastaReader.readFasta(fnaFile);
                }
            }
        }
        Map<File, TaxTree.TaxIdNode> additionalMap = additionalGoal.get();
        for (File additionalFasta : additionalMap.keySet()) {
            TaxTree.TaxIdNode node = additionalMap.get(additionalFasta);
            if (taxNodesGoal.get().contains(node)) {
                fastaReader.ignoreAccessionMap(node);
                fastaReader.readFasta(additionalFasta);
            }
        }
    }

    @Override
    public boolean isWeakDependency(Goal<GSProject> toGoal) {
        if (toGoal == fnaFilesGoal) {
            return true;
        }
        return super.isWeakDependency(toGoal);
    }
}
