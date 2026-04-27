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

import it.unimi.dsi.fastutil.objects.Object2LongMap;
import org.metagene.genestrip.GSGoalKey;
import org.metagene.genestrip.GSMaker;
import org.metagene.genestrip.finertree.cluster.DendrogramNode;
import org.metagene.genestrip.finertree.goals.*;
import org.metagene.genestrip.goals.*;
import org.metagene.genestrip.goals.refseq.RefSeqFnaFilesDownloadGoal;
import org.metagene.genestrip.goals.refseq.SVGTaxTreeGoal;
import org.metagene.genestrip.goals.refseq.StoreDBGoal;
import org.metagene.genestrip.io.StreamingResourceStream;
import org.metagene.genestrip.make.FileGoal;
import org.metagene.genestrip.make.FileListGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.match.MatchingResult;
import org.metagene.genestrip.refseq.AccessionMap;
import org.metagene.genestrip.refseq.RefSeqCategory;
import org.metagene.genestrip.store.Database;
import org.metagene.genestrip.tax.SmallTaxTree;
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
        KMerIndexBloomGoal<P> bloomGoal = new KMerIndexBloomGoal(project, getExecutionContext(project),
                categoriesGoal, taxNodesGoal, taxTreeGoal, fnaFilesGoal, additionalGoal, accessionMapGoal, storeGoal);
        registerGoal(bloomGoal);

        StoreKMerIndexGoal<P> storeKMerIndexGoal = new StoreKMerIndexGoal(project, bloomGoal);
        registerGoal(storeKMerIndexGoal);

        LoadKMerIndexGoal<P> loadKMerIndexGoal = new LoadKMerIndexGoal(project, bloomGoal, storeKMerIndexGoal);
        registerGoal(loadKMerIndexGoal);

        KMerIntersectCountGoal<P> intersectCountGoal = new KMerIntersectCountGoal(project, storeGoal, loadKMerIndexGoal);
        registerGoal(intersectCountGoal);

        KMerIntersectCSVGoal<P> csvGoal = new KMerIntersectCSVGoal(project, storeGoal, intersectCountGoal);
        registerGoal(csvGoal);

        ObjectGoal<Map<SmallTaxTree.SmallTaxIdNode, DendrogramNode>, P> dendrogramGoal = new DendrogramGoal(project, intersectCountGoal);
        registerGoal(dendrogramGoal);

        DengrogramLaTeXGoal<P> laTeXGoal = new DengrogramLaTeXGoal(project, storeGoal, dendrogramGoal, projectSetupGoal);
        registerGoal(laTeXGoal);

        UpdateStoreGoal<P> updateStoreGoal = new UpdateStoreGoal(project, storeGoal, dendrogramGoal, loadKMerIndexGoal);
        registerGoal(updateStoreGoal);

        StoreDBGoal<P> storeUpdatedDBGoal = new StoreDBGoal(project, FTGoalKey.FTDB,
                project.getOutputFile(FTGoalKey.FTDB.getName(), P.GSFileType.DB, false), updateStoreGoal);
        registerGoal(storeUpdatedDBGoal);

        LoadDBGoal<P> loadFTDBGoal = new LoadDBGoal(project, FTGoalKey.LOAD_FTDB, updateStoreGoal, storeUpdatedDBGoal);
        registerGoal(loadFTDBGoal);

        FTDBInfoGoal infoGoal = new FTDBInfoGoal(project, loadFTDBGoal);
        registerGoal(infoGoal);

        FileGoal<P> allInOneLaTeXGoal = new AllInOneLaTeXGoal(project, laTeXGoal, projectSetupGoal);
        registerGoal(allInOneLaTeXGoal);

        ObjectGoal<Map<String, StreamingResourceStream>, P> fastqMapTransfGoal = (ObjectGoal<Map<String, StreamingResourceStream>, P>) getGoal(GSGoalKey.FASTQ_MAP_TRANSFORM);
        FastqDownloadsGoal<P> fastqDownloadsGoal = (FastqDownloadsGoal<P>) getGoal(GSGoalKey.FASTQ_DOWNLOAD);

        ObjectGoal<Map<String, MatchingResult>, P> ftmatchResGoal = new MatchResultGoal(getProject(), FTGoalKey.FTMATCHRES, fastqMapTransfGoal, loadFTDBGoal,
                getExecutionContext(getProject()), projectSetupGoal, fastqDownloadsGoal);
        registerGoal(ftmatchResGoal);

        Goal<P> ftmatchGoal = new MatchGoal(project, FTGoalKey.FTMATCH, fastqMapTransfGoal, ftmatchResGoal, projectSetupGoal);
        registerGoal(ftmatchGoal);

        ObjectGoal<Set<SmallTaxTree.SmallTaxIdNode>, P> db2fastqTaxNodesGoal = (ObjectGoal<Set<SmallTaxTree.SmallTaxIdNode>, P>) getGoal(GSGoalKey.DB2FASTQ_TAXIDS);
        Goal<P> db2fastqGoal = new DB2FastqGoal(project, FTGoalKey.FTDB2FASTQ, db2fastqTaxNodesGoal, loadFTDBGoal, projectSetupGoal);
        registerGoal(db2fastqGoal);

        SVGTaxTreeGoal<P> svgTaxTreeGoal = new SVGTaxTreeGoal<P>(project, FTGoalKey.FT_SVG_TAX_TREE, loadFTDBGoal, projectSetupGoal);
        registerGoal(svgTaxTreeGoal);

        ObjectGoal<Object2LongMap<String>, P> ftKmersPerTaxGoal = new KMersPerDBTaxidGoal<>(project, FTGoalKey.FT_KMERS_PER_TAXID, getExecutionContext(project),
                categoriesGoal, taxNodesGoal, fnaFilesGoal, additionalGoal, accessionMapGoal, loadFTDBGoal, taxTreeGoal /* taxTreeGoal is only REQUIRED so that the tree is not dropped to early! */);
        registerGoal(ftKmersPerTaxGoal);

        Goal<P> ftConsGoal = new FTConsistencyGoal<>(project, FTGoalKey.FT_CONSISTENCY, loadFTDBGoal, ftKmersPerTaxGoal);
        registerGoal(ftConsGoal);

        ObjectGoal<Object2LongMap<String>, P> kmersPerTaxGoal = new KMersPerDBTaxidGoal<>(project, FTGoalKey.KMERS_PER_TAXID, getExecutionContext(project),
                categoriesGoal, taxNodesGoal, fnaFilesGoal, additionalGoal, accessionMapGoal, storeGoal, taxTreeGoal /* taxTreeGoal is only REQUIRED so that the tree is not dropped to early! */);
        registerGoal(kmersPerTaxGoal);

        Goal<P> dbConsGoal = new FTConsistencyGoal<>(project, FTGoalKey.DB_CONSISTENCY, storeGoal, kmersPerTaxGoal);
        registerGoal(dbConsGoal);

        Goal<P> clearGoal = getGoal(GSGoalKey.CLEAR);
        Goal<P> ftclearGoal = new FileListGoal<P>(project, FTGoalKey.FTCLEAR, Arrays
                .asList(project.getTeXDir()), clearGoal) {
            @Override
            public boolean isMade() {
                return false;
            }

            @Override
            protected void makeFile(File file) throws IOException {
            }

            @Override
            protected void doMakeThis() {
                doCleanThis();
            }
        };
        registerGoal(ftclearGoal);
    }

    @Override
    protected LoadDBGoal<P> getLoadDBGoal() {
        if (isUseFTDBForAPI()) {
            return (LoadDBGoal) getGoal(FTGoalKey.LOAD_FTDB);
        }
        else {
            return super.getLoadDBGoal();
        }
    }
}
