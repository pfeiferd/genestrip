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
package org.metagene.genestrip.finertree;

import org.metagene.genestrip.GSGoalKey;
import org.metagene.genestrip.GSMaker;
import org.metagene.genestrip.finertree.goals.*;
import org.metagene.genestrip.goals.refseq.RefSeqFnaFilesDownloadGoal;
import org.metagene.genestrip.make.FileListGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.refseq.AccessionMap;
import org.metagene.genestrip.refseq.RefSeqCategory;
import org.metagene.genestrip.store.Database;
import org.metagene.genestrip.tax.TaxTree;

import java.io.File;
import java.io.IOException;
import java.util.*;

public class FinerTreeMaker<P extends FTProject> extends GSMaker<P> {
    private boolean useFTDBForAPI;

    public FinerTreeMaker(P project) {
        super(project);
        setUseFTDBForAPI(false);
    }

    public void setUseFTDBForAPI(boolean useFTDBForAPI) {
        this.useFTDBForAPI = useFTDBForAPI;
    }

    public boolean isUseFTDBForAPI() {
        return useFTDBForAPI;
    }

    @Override
    protected void createGoals() {
        super.createGoals();

        P project = getProject();

        List<File> projectDirs = Arrays.asList(project.getTeXDir());
        Goal<P> projectSetupGoal = new FileListGoal<P>(project, FTGoalKey.FTSETUP, projectDirs,
                getGoal(GSGoalKey.SETUP)) {
            @Override
            protected void makeFile(File file) throws IOException {
                file.mkdir();
            }

            @Override
            public boolean isAllowTransitiveClean() {
                return false;
            }
        };
        registerGoal(projectSetupGoal);

        ObjectGoal<Set<RefSeqCategory>, P> categoriesGoal = (ObjectGoal<Set<RefSeqCategory>, P>) getGoal(GSGoalKey.CATEGORIES);
        ObjectGoal<Set<TaxTree.TaxIdNode>, P> taxNodesGoal = (ObjectGoal<Set<TaxTree.TaxIdNode>, P>) getGoal(GSGoalKey.TAXNODES);
        ObjectGoal<TaxTree, P> taxTreeGoal = (ObjectGoal<TaxTree, P>) getGoal(GSGoalKey.TAXTREE);
        RefSeqFnaFilesDownloadGoal fnaFilesGoal = (RefSeqFnaFilesDownloadGoal) getGoal(GSGoalKey.REFSEQFNA);
        ObjectGoal<Map<File, TaxTree.TaxIdNode>, P> additionalGoal = (ObjectGoal<Map<File, TaxTree.TaxIdNode>, P>) getGoal(GSGoalKey.ADD_FASTAS);
        ObjectGoal<AccessionMap, P> accessionMapGoal = (ObjectGoal<AccessionMap, P>) getGoal(GSGoalKey.ACCMAP);
        ObjectGoal<Database, P> storeGoal = (ObjectGoal<Database, P>) getGoal(GSGoalKey.LOAD_DB);

        ObjectGoal<Map<String, DBQualityCountsGoal.Counts>, P> kmersPerTaxGoal = new DBQualityCountsGoal<>(project, FTGoalKey.DB_QUALITY_COUNTS, getExecutionContext(project),
                categoriesGoal, taxNodesGoal, fnaFilesGoal, additionalGoal, accessionMapGoal, storeGoal, taxTreeGoal /* taxTreeGoal is only REQUIRED so that the tree is not dropped too early! */);
        registerGoal(kmersPerTaxGoal);

        Goal<P> dbQualityGoal = new DBQualityCSVGoal<>(project, FTGoalKey.DB_QUALITY, storeGoal, kmersPerTaxGoal);
        registerGoal(dbQualityGoal);
    }
}
