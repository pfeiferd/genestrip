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
import org.metagene.genestrip.tax.TaxNodeSelection;

import java.io.File;
import java.io.IOException;
import java.util.*;

/**
 * {@link GSMaker} for the finer-tree (FT) extension. Builds the standard Genestrip goal graph and adds
 * the FT-specific setup and database-quality goals.
 *
 * @param <P> the concrete FT project type
 */
public class FinerTreeMaker<P extends FTProject> extends GSMaker<P> {
    private boolean useFTDBForAPI;

    /**
     * Creates a finer-tree maker for the given project.
     *
     * @param project the project the goals are built for
     */
    public FinerTreeMaker(P project) {
        super(project);
        setUseFTDBForAPI(false);
    }

    /**
     * Sets whether the FT database should be used for API access.
     *
     * @param useFTDBForAPI whether the FT database is used for API access
     */
    public void setUseFTDBForAPI(boolean useFTDBForAPI) {
        this.useFTDBForAPI = useFTDBForAPI;
    }

    /**
     * Returns whether the FT database is used for API access.
     *
     * @return whether the FT database is used for API access
     */
    public boolean isUseFTDBForAPI() {
        return useFTDBForAPI;
    }

    /**
     * Registers the standard Genestrip goals and additionally the FT setup goal and the
     * database-quality counting and CSV goals.
     */
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
        ObjectGoal<TaxNodeSelection, P> taxNodesGoal = (ObjectGoal<TaxNodeSelection, P>) getGoal(GSGoalKey.TAXNODES);
        ObjectGoal<TaxTree, P> taxTreeGoal = (ObjectGoal<TaxTree, P>) getGoal(GSGoalKey.TAXTREE);
        RefSeqFnaFilesDownloadGoal fnaFilesGoal = (RefSeqFnaFilesDownloadGoal) getGoal(GSGoalKey.REFSEQFNA);
        ObjectGoal<Map<File, TaxTree.TaxIdNode>, P> additionalGoal = (ObjectGoal<Map<File, TaxTree.TaxIdNode>, P>) getGoal(GSGoalKey.ADD_FASTAS);
        ObjectGoal<AccessionMap, P> accessionMapGoal = (ObjectGoal<AccessionMap, P>) getGoal(GSGoalKey.ACCMAP);
        ObjectGoal<Database, P> storeGoal = (ObjectGoal<Database, P>) getGoal(GSGoalKey.LOAD_DB);
        // The sizing goal is registered like any other but depended upon lazily: being an ObjectGoal it
        // is never made on its own account, only when the bloom goal asks it for its value, which it
        // does exactly when `kmerIndexSizing' calls for the estimate. Configured to the conservative
        // bound, this goal never runs and nothing is read twice.
        KMerIndexSizeGoal<P> indexSizeGoal = new KMerIndexSizeGoal(project, getExecutionContext(project),
                categoriesGoal, taxNodesGoal, taxTreeGoal, fnaFilesGoal, additionalGoal, accessionMapGoal, storeGoal);
        registerGoal(indexSizeGoal);

        KMerIndexBloomGoal<P> bloomGoal = new KMerIndexBloomGoal(project, getExecutionContext(project),
                categoriesGoal, taxNodesGoal, taxTreeGoal, fnaFilesGoal, additionalGoal, accessionMapGoal, storeGoal,
                indexSizeGoal);
        registerGoal(bloomGoal);

        StoreKMerIndexGoal<P> storeKMerIndexGoal = new StoreKMerIndexGoal(project, bloomGoal);
        registerGoal(storeKMerIndexGoal);

        LoadKMerIndexGoal<P> loadKMerIndexGoal = new LoadKMerIndexGoal(project, bloomGoal, storeKMerIndexGoal);
        registerGoal(loadKMerIndexGoal);

        KMerIntersectCountGoal<P> intersectCountGoal = new KMerIntersectCountGoal(project, storeGoal, loadKMerIndexGoal);
        registerGoal(intersectCountGoal);

        KMerIntersectCSVGoal<P> csvGoal = new KMerIntersectCSVGoal(project, storeGoal, intersectCountGoal);
        registerGoal(csvGoal);

        KMerBranchHistoGoal<P> branchHistoGoal = new KMerBranchHistoGoal(project, storeGoal, loadKMerIndexGoal);
        registerGoal(branchHistoGoal);

        KMerBranchHistoCSVGoal<P> branchHistoCSVGoal = new KMerBranchHistoCSVGoal(project, storeGoal, branchHistoGoal);
        registerGoal(branchHistoCSVGoal);

        KMerBranchHistoRankCSVGoal<P> branchHistoRankCSVGoal = new KMerBranchHistoRankCSVGoal(project, storeGoal, branchHistoGoal);
        registerGoal(branchHistoRankCSVGoal);

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

        ObjectGoal<Long, P> ftQualitySizeGoal = new DBQualitySizeGoal<>(project, FTGoalKey.FT_QUALITY_SIZE, getExecutionContext(project),
                categoriesGoal, taxNodesGoal, fnaFilesGoal, additionalGoal, accessionMapGoal, loadFTDBGoal, taxTreeGoal);
        registerGoal(ftQualitySizeGoal);

        ObjectGoal<Map<String, DBQualityCountsGoal.Counts>, P> ftKmersPerTaxGoal = new DBQualityCountsGoal<>(project, FTGoalKey.FT_QUALITY_COUNTS, getExecutionContext(project),
                categoriesGoal, taxNodesGoal, fnaFilesGoal, additionalGoal, accessionMapGoal, loadFTDBGoal, ftQualitySizeGoal, taxTreeGoal /* taxTreeGoal is only REQUIRED so that the tree is not dropped too early! */);
        registerGoal(ftKmersPerTaxGoal);

        Goal<P> ftQualityGoal = new DBQualityCSVGoal<>(project, FTGoalKey.FT_QUALITY, loadFTDBGoal, ftKmersPerTaxGoal);
        registerGoal(ftQualityGoal);

        ObjectGoal<Long, P> dbQualitySizeGoal = new DBQualitySizeGoal<>(project, FTGoalKey.DB_QUALITY_SIZE, getExecutionContext(project),
                categoriesGoal, taxNodesGoal, fnaFilesGoal, additionalGoal, accessionMapGoal, storeGoal, taxTreeGoal);
        registerGoal(dbQualitySizeGoal);

        ObjectGoal<Map<String, DBQualityCountsGoal.Counts>, P> kmersPerTaxGoal = new DBQualityCountsGoal<>(project, FTGoalKey.DB_QUALITY_COUNTS, getExecutionContext(project),
                categoriesGoal, taxNodesGoal, fnaFilesGoal, additionalGoal, accessionMapGoal, storeGoal, dbQualitySizeGoal, taxTreeGoal /* taxTreeGoal is only REQUIRED so that the tree is not dropped too early! */);
        registerGoal(kmersPerTaxGoal);

        Goal<P> dbQualityGoal = new DBQualityCSVGoal<>(project, FTGoalKey.DB_QUALITY, storeGoal, kmersPerTaxGoal);
        registerGoal(dbQualityGoal);

        Goal<P> ftAll = new Goal<P>(project, FTGoalKey.FTGENALL, getGoal(GSGoalKey.GENALL), infoGoal) {
            @Override
            public boolean isMade() {
                return false;
            }

            @Override
            protected void doMakeThis() {
            }
        };
        registerGoal(ftAll);
        // Supersedes GSGoalKey.GENALL, which GSMaker registers as its default goal.
        setDefaultGoal(ftAll);

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

    /**
     * Returns the goal used to load the database for API access, choosing the FT database
     * ({@link FTGoalKey#LOAD_FTDB}) when {@link #isUseFTDBForAPI()} is set and otherwise the standard
     * Genestrip load goal.
     *
     * @return the load-database goal used for API access
     */
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
