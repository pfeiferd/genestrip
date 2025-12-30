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
package org.metagene.genestrip;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.metagene.genestrip.GSProject.GSFileType;
import org.metagene.genestrip.goals.*;
import org.metagene.genestrip.goals.genbank.AssemblyFileDownloadGoal;
import org.metagene.genestrip.goals.genbank.FastaFilesFromGenbankGoal;
import org.metagene.genestrip.goals.genbank.FastaFilesGenbankDownloadGoal;
import org.metagene.genestrip.goals.genbank.TaxNodesFromGenbankGoal;
import org.metagene.genestrip.goals.kraken.KrakenResCountGoal;
import org.metagene.genestrip.goals.kraken.KrakenResFileGoal;
import org.metagene.genestrip.goals.refseq.*;
import org.metagene.genestrip.io.StreamingResourceStream;
import org.metagene.genestrip.make.FileGoal;
import org.metagene.genestrip.make.FileListGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.Maker;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.match.MatchingResult;
import org.metagene.genestrip.refseq.AccessionMap;
import org.metagene.genestrip.refseq.RefSeqCategory;
import org.metagene.genestrip.tax.SmallTaxTree.SmallTaxIdNode;
import org.metagene.genestrip.tax.TaxTree;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;

public class GSMaker<P extends GSProject> extends Maker<P> {
    private ExecutionContext executionContext;

    public GSMaker(P project) {
        super(project);
    }

    public void dumpAll() {
        dump();
        if (executionContext != null) {
            executionContext.dump();
        }
    }

    protected ExecutionContext createExecutionContext(Thread mainThread, P project) {
        return new DefaultExecutionContext(mainThread, project.intConfigValue(GSConfigKey.THREADS),
                project.longConfigValue(GSConfigKey.LOG_PROGRESS_UPDATE_CYCLE));
    }

    protected ExecutionContext getExecutionContext(P project) {
        if (executionContext == null) {
            executionContext = createExecutionContext(getMainThread(), project);
        }
        return executionContext;
    }

    protected Thread getMainThread() {
        return Thread.currentThread();
    }

    protected void createGoals() {
        P project = getProject();
        // Show goals...

        Goal<P> showAllGoals = new Goal<P>(project, GSGoalKey.SHOWALL) {
            @Override
            public boolean isMade() {
                return false;
            }

            @Override
            protected void doMakeThis() {
                List<String> res = new ArrayList<>();
                for (GSGoalKey key : GSGoalKey.values()) {
                    if (!(getGoal(key) instanceof ObjectGoal<?, ?>)) {
                        res.add(key.getName());
                    }
                }
                System.out.println(res);
            }
        };
        registerGoal(showAllGoals);

        Goal<P> showGoals = new Goal<P>(project, GSGoalKey.SHOW) {
            @Override
            public boolean isMade() {
                return false;
            }

            @Override
            protected void doMakeThis() {
                List<String> res = new ArrayList<>();
                for (GSGoalKey key : GSGoalKey.values()) {
                    if (key.isForUser()) {
                        res.add(key.getName());
                    }
                }
                System.out.println(res);
            }
        };
        registerGoal(showGoals);

        // Common setup

        Goal<P> commonSetupGoal = new FileListGoal<P>(project, GSGoalKey.COMMON_SETUP,
                Arrays.asList(project.getCommon().getCommonDir(), project.getCommon().getRefSeqDir(),
                        project.getCommon().getGenbankDir(), project.getCommon().getFastqDir(),
                        project.getCommon().getFastaDir())) {
            @Override
            protected void makeFile(File file) throws IOException {
                file.mkdir();
            }

            @Override
            protected List<File> getFilesToClean() {
                List<File> files = new ArrayList<>(getFiles());
                files.remove(project.getCommon().getFastqDir());
                files.remove(project.getCommon().getFastaDir());
                return files;
            }

            @Override
            public boolean isAllowTransitiveClean() {
                return false;
            }
        };
        registerGoal(commonSetupGoal);

        // Taxonomy

        Goal<P> taxDBGoal = new TaxIdFileDownloadGoal(project, commonSetupGoal);
        registerGoal(taxDBGoal);

        ObjectGoal<TaxTree, P> taxTreeGoal = new ObjectGoal<TaxTree, P>(project, GSGoalKey.TAXTREE,
                taxDBGoal) {
            @Override
            protected void doMakeThis() {
                set(new TaxTree(getProject().getCommon().getCommonDir(), booleanConfigValue(GSConfigKey.PROGRESS_BAR)));
            }
        };
        registerGoal(taxTreeGoal);

        ObjectGoal<Set<TaxIdNode>, P> taxNodesGoal = new TaxNodesGoal(project, taxTreeGoal);
        registerGoal(taxNodesGoal);

        // General stuff from RefSeq

        FileGoal<P> releaseNumberGoal = new RefSeqRNumDownloadGoal(project, commonSetupGoal);
        registerGoal(releaseNumberGoal);

        RefSeqCatalogDownloadGoal refSeqCatalogGoal = new RefSeqCatalogDownloadGoal(project, releaseNumberGoal,
                commonSetupGoal);
        registerGoal(refSeqCatalogGoal);

        CheckSumMapGoal checkSumMapGoal = new CheckSumMapGoal(project, refSeqCatalogGoal);
        registerGoal(checkSumMapGoal);

        // Create or clear project directories
        List<File> projectDirs = Arrays.asList(project.getFastaDir(), project.getFastqDir(), project.getDBDir(),
                project.getKrakenOutDir(), project.getResultsDir(), project.getLogDir());

        Goal<P> projectSetupGoal = new FileListGoal<P>(project, GSGoalKey.SETUP, projectDirs,
                commonSetupGoal) {
            @Override
            protected List<File> getFilesToClean() {
                List<File> files = new ArrayList<>(getFiles());
                files.remove(project.getFastaDir());
                files.remove(project.getFastqDir());
                return files;
            }

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

        Goal<P> clearGoal = new FileListGoal<P>(project, GSGoalKey.CLEAR, Arrays
                .asList(project.getDBDir(), project.getKrakenOutDir(), project.getResultsDir(), project.getLogDir())) {
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
        registerGoal(clearGoal);

        // Download genomic data for project

        ObjectGoal<Set<RefSeqCategory>, P> categoriesGoal = new CategoriesGoal(project, projectSetupGoal);
        registerGoal(categoriesGoal);

        ObjectGoal<CheckRefSeqRNumGoal.Result, P> checkRefSeqRNumGoal = new CheckRefSeqRNumGoal(project, releaseNumberGoal);
        registerGoal(checkRefSeqRNumGoal);

        RefSeqFnaFilesDownloadGoal refSeqFnaFilesGoal = new RefSeqFnaFilesDownloadGoal(project, categoriesGoal,
                refSeqCatalogGoal, checkSumMapGoal, checkRefSeqRNumGoal);
        registerGoal(refSeqFnaFilesGoal);

        ObjectGoal<Integer, P> accessMapSizeGoal = new AccessionMapSizeGoal(project, categoriesGoal,
                refSeqCatalogGoal);
        registerGoal(accessMapSizeGoal);

        ObjectGoal<AccessionMap, P> accessionMapGoal = new AccessionMapGoal(project, categoriesGoal, taxTreeGoal,
                refSeqCatalogGoal, accessMapSizeGoal);
        registerGoal(accessionMapGoal);

        // Genbank related additional fastas:
        ObjectGoal<Set<TaxIdNode>, P> taxNodesFromGenBank = new TaxNodesFromGenbankGoal(project, categoriesGoal,
                taxNodesGoal, accessionMapGoal);
        registerGoal(taxNodesFromGenBank);

        FileGoal<P> assemblyFileDownloadGoal = new AssemblyFileDownloadGoal(project, commonSetupGoal);
        registerGoal(assemblyFileDownloadGoal);

        FastaFilesFromGenbankGoal<P> fastaFilesFromGenbankGoal = new FastaFilesFromGenbankGoal(project, taxTreeGoal,
                assemblyFileDownloadGoal, taxNodesFromGenBank);
        registerGoal(fastaFilesFromGenbankGoal);

        FastaFilesGenbankDownloadGoal<P> fastaFilesGenbankDownloadGoal = new FastaFilesGenbankDownloadGoal(project,
                fastaFilesFromGenbankGoal);
        registerGoal(fastaFilesGenbankDownloadGoal);
        // end of Genbank stuff.

        AdditionalDownloadsGoal<P> additionalDownloadsGoal = new AdditionalDownloadsGoal(project, projectSetupGoal);
        registerGoal(additionalDownloadsGoal);

        ObjectGoal<Map<File, TaxIdNode>, P> additionalFastasGoal = new AdditionalFastasGoal(project,
                taxTreeGoal, fastaFilesFromGenbankGoal, fastaFilesGenbankDownloadGoal, projectSetupGoal,
                additionalDownloadsGoal);
        registerGoal(additionalFastasGoal);

        // Create database and bloom filter

        FillSizeGoal<P> fillSizeGoal = new FillSizeGoal(project, getExecutionContext(project), categoriesGoal, taxNodesGoal, refSeqFnaFilesGoal,
                additionalFastasGoal, accessionMapGoal);
        registerGoal(fillSizeGoal);

        ObjectGoal<Long, P> fillBloomGoal = new FillBloomFilterGoal(project, getExecutionContext(project), categoriesGoal,
                taxNodesGoal, refSeqFnaFilesGoal, additionalFastasGoal, accessionMapGoal, fillSizeGoal);
        registerGoal(fillBloomGoal);

        FillDBGoal<P> fillDBGoal = new FillDBGoal(project, getExecutionContext(project), categoriesGoal, taxNodesGoal, taxTreeGoal, refSeqFnaFilesGoal,
                additionalFastasGoal, accessionMapGoal, fillBloomGoal, projectSetupGoal);
        registerGoal(fillDBGoal);

        StoreDBGoal<P> storeTempDBGoal = new StoreDBGoal<P>(project, GSGoalKey.TEMPDB,
                project.getOutputFile(GSGoalKey.TEMPDB.getName(), GSFileType.DB, false), fillDBGoal, projectSetupGoal) {
            @Override
            protected void dependentMade(Goal<P> goal) {
                super.dependentMade(goal);
                if (booleanConfigValue(GSConfigKey.REMOVE_TEMP_DB)) {
                    // Remove temp file if not required any more...
                    if (GSGoalKey.DB.equals(goal.getKey())) {
                        cleanThis();
                    }
                }
            }
        };
        registerGoal(storeTempDBGoal);

        FilledDBGoal<P> filledDBGoal = new FilledDBGoal(project, fillDBGoal, storeTempDBGoal);
        registerGoal(filledDBGoal);

        Goal<P> tempDbInfoGoal = new DBInfoGoal<P>(true, project, filledDBGoal, projectSetupGoal) {
            @Override
            protected void dependentMade(Goal<P> goal) {
                super.dependentMade(goal);
                if (booleanConfigValue(GSConfigKey.REMOVE_TEMP_DB)) {
                    // Remove temp file if not required any more...
                    if (GSGoalKey.DB.equals(goal.getKey())) {
                        cleanThis();
                    }
                }
            }
        };
        registerGoal(tempDbInfoGoal);

        DBGoal<P> updateDBGoal = new DBGoal(project, getExecutionContext(project), categoriesGoal, taxNodesGoal, taxTreeGoal,
                refSeqFnaFilesGoal, additionalFastasGoal, accessionMapGoal, filledDBGoal, tempDbInfoGoal, projectSetupGoal);
        registerGoal(updateDBGoal);

        // We weed the storeTempDBGoal and tempDbInfoGoal as an explicit dependency here so that associated intermediate files get
        // automatically deleted when final database got stored.
        StoreDBGoal<P> storeDBGoal = new StoreDBGoal(project, GSGoalKey.DB, project.getDBFile(), updateDBGoal,
                projectSetupGoal,
                // Keep these two dependencies:
                storeTempDBGoal, tempDbInfoGoal);
        registerGoal(storeDBGoal);

        LoadDBGoal<P> loadDBGoal = new LoadDBGoal(project, GSGoalKey.LOAD_DB, updateDBGoal, storeDBGoal);
        registerGoal(loadDBGoal);

        BloomIndexGoal<P> bloomIndexGoal = new BloomIndexGoal(project, loadDBGoal, projectSetupGoal);
        registerGoal(bloomIndexGoal);

        StoreIndexGoal<P> storeIndexGoal = new StoreIndexGoal(project, bloomIndexGoal, projectSetupGoal);
        registerGoal(storeIndexGoal);

        LoadIndexGoal<P> bloomIndexedGoal = new LoadIndexGoal(project, bloomIndexGoal, storeIndexGoal, projectSetupGoal);
        registerGoal(bloomIndexedGoal);

        // Manage the database

        Goal<P> dbInfoGoal = new DBInfoGoal(project, loadDBGoal, projectSetupGoal);
        registerGoal(dbInfoGoal);

        Goal<P> all = new Goal<P>(project, GSGoalKey.GENALL, dbInfoGoal, storeIndexGoal) {
            @Override
            public boolean isMade() {
                return false;
            }

            @Override
            protected void doMakeThis() {
            }
        };
        registerDefaultGoal(all);

        ObjectGoal<Set<SmallTaxIdNode>, P> db2fastqTaxNodesGoal = new DB2FastqTaxNodesGoal(project, loadDBGoal,
                projectSetupGoal);
        registerGoal(db2fastqTaxNodesGoal);

        Goal<P> db2fastqGoal = new DB2FastqGoal(project, GSGoalKey.DB2FASTQ, db2fastqTaxNodesGoal, loadDBGoal, projectSetupGoal);
        registerGoal(db2fastqGoal);

        // Use database and bloom filter

        ObjectGoal<Map<String, StreamingResourceStream>, P> fastqMapGoal = new FastqMapGoal(project, true,
                projectSetupGoal);
        registerGoal(fastqMapGoal);

        ObjectGoal<Map<String, StreamingResourceStream>, P> fastqMapTransfGoal = new FastqMapTransformGoal(
                project, true, fastqMapGoal, projectSetupGoal);
        registerGoal(fastqMapTransfGoal);

        Goal<P> extractGoal = new ExtractGoal(project, fastqMapTransfGoal, getExecutionContext(project));
        registerGoal(extractGoal);

        FastqDownloadsGoal<P> fastqDownloadsGoal = new FastqDownloadsGoal(project, true, fastqMapGoal, fastqMapTransfGoal,
                projectSetupGoal);
        registerGoal(fastqDownloadsGoal);

        Goal<P> filterGoal = new FilterGoal(project, fastqMapTransfGoal, bloomIndexedGoal,
                getExecutionContext(project), projectSetupGoal, fastqDownloadsGoal);
        registerGoal(filterGoal);

        ObjectGoal<Map<String, MatchingResult>, P> matchResGoal = new MatchResultGoal(getProject(), GSGoalKey.MATCHRES, fastqMapTransfGoal, loadDBGoal,
                getExecutionContext(getProject()), projectSetupGoal, fastqDownloadsGoal);
        registerGoal(matchResGoal);

        Goal<P> matchGoal = new MatchGoal(project, GSGoalKey.MATCH, fastqMapTransfGoal, matchResGoal, projectSetupGoal);
        registerGoal(matchGoal);

        ObjectGoal<Map<String, MatchingResult>, P> matchReslrGoal = new MatchResultGoal(getProject(), GSGoalKey.MATCHRESLR, fastqMapTransfGoal, loadDBGoal,
                getExecutionContext(getProject()), projectSetupGoal, fastqDownloadsGoal);
        registerGoal(matchReslrGoal);

        Goal<P> matchlrGoal = new MatchGoal(project, GSGoalKey.MATCHLR, fastqMapTransfGoal, matchReslrGoal, projectSetupGoal);
        registerGoal(matchlrGoal);

        ObjectGoal<Map<String, StreamingResourceStream>, P> fastaMapGoal = new FastqMapGoal(project, false,
                projectSetupGoal);
        registerGoal(fastaMapGoal);

        // Fasta to fastq

        ObjectGoal<Map<String, StreamingResourceStream>, P> fastaMapTransfGoal = new FastqMapTransformGoal(
                project, false, fastaMapGoal, projectSetupGoal);
        registerGoal(fastaMapTransfGoal);

        FastqDownloadsGoal<P> fastaDownloadsGoal = new FastqDownloadsGoal(project, false, fastaMapGoal, fastaMapTransfGoal,
                projectSetupGoal);
        registerGoal(fastaDownloadsGoal);

        Goal<P> fastq2fastaGoal = new Fasta2FastqGoal(project, GSGoalKey.FASTA2FASTQ, fastaMapTransfGoal, projectSetupGoal, fastaDownloadsGoal);
        registerGoal(fastq2fastaGoal);

        // Use kraken
        KrakenResCountGoal<P> krakenResCountGoal = new KrakenResCountGoal(project, fastqMapTransfGoal, taxNodesGoal,
                projectSetupGoal, fastqDownloadsGoal);
        registerGoal(krakenResCountGoal);

        KrakenResFileGoal krakenResFileGoal = new KrakenResFileGoal(project, fastqMapTransfGoal, taxNodesGoal,
                krakenResCountGoal, projectSetupGoal, fastqDownloadsGoal);
        registerGoal(krakenResFileGoal);
    }

    protected void createGoalDependingOnLoadDB() {

    }

    public MatchingResult cleanMatch(boolean lr, String key, String... pathsOrURLs) {
        return match(lr, true, key, pathsOrURLs);
    }

    public MatchingResult match(boolean lr, String key, String... pathsOrURLs) {
        return match(lr, false, key, pathsOrURLs);
    }

    public MatchingResult match(boolean lr, boolean clean, String key, String... pathsOrURLs) {
        MatchResultGoal<P> matchResGoal = createGoalChainForMatchResult(lr, key, pathsOrURLs);
        MatchGoal<P> matchGoal = new MatchGoal(matchResGoal.getProject(), (lr ? GSGoalKey.MATCHLR : GSGoalKey.MATCH), matchResGoal.getFastqMapGoal(), matchResGoal);
        if (clean) {
            matchGoal.cleanThis();
        }
        // Store it here as it will be cleared via "allDependentsMade()" in matchGoal.make()...
        MatchingResult res = matchResGoal.get().get(key);
        matchGoal.make();
        return res;
    }

    public MatchingResult matchResult(boolean lr, String key, String... pathsOrURLs) {
        return createGoalChainForMatchResult(lr, key, pathsOrURLs).get().get(key);
    }

    protected MatchResultGoal<P> createGoalChainForMatchResult(boolean lr, String key, String... pathsOrURLs) {
        ObjectGoal<Map<String, StreamingResourceStream>, P> fastqMapGoal = new FastqMapGoal(getProject(), true,
                getGoal(GSGoalKey.SETUP)) {
            @Override
            protected void doMakeThis() {
                Map<String, StreamingResourceStream> map = createFastqMap(key, pathsOrURLs, null, null, null);
                set(map);
                if (getLogger().isInfoEnabled()) {
                    getLogger().info("Derived fastq map: " + map);
                }
            }
        };

        ObjectGoal<Map<String, StreamingResourceStream>, P> fastqMapTransfGoal = new FastqMapTransformGoal(
                getProject(), true, fastqMapGoal, getGoal(GSGoalKey.SETUP));

        FastqDownloadsGoal<P> fastqDownloadsGoal = new FastqDownloadsGoal(getProject(), true, fastqMapGoal, fastqMapTransfGoal,
                getGoal(GSGoalKey.SETUP));

        LoadDBGoal<P> loadDBGoal = getLoadDBGoal();

        return new MatchResultGoal<P>(getProject(), (lr ? GSGoalKey.MATCHRESLR : GSGoalKey.MATCHRES), fastqMapTransfGoal, loadDBGoal,
                getExecutionContext(getProject()), getGoal(GSGoalKey.SETUP), fastqDownloadsGoal);
    }

    protected LoadDBGoal<P> getLoadDBGoal() {
        return (LoadDBGoal) getGoal(GSGoalKey.LOAD_DB);
    }

    protected FilterGoal<P> createGoalChainForFilter(String key, String... pathsOrURLs) {
        ObjectGoal<Map<String, StreamingResourceStream>, P> fastqMapGoal = new FastqMapGoal(getProject(), true,
                getGoal(GSGoalKey.SETUP)) {
            @Override
            protected void doMakeThis() {
                Map<String, StreamingResourceStream> map = createFastqMap(key, pathsOrURLs, null, null, null);
                set(map);
                if (getLogger().isInfoEnabled()) {
                    getLogger().info("Derived fastq map: " + map);
                }
            }
        };

        ObjectGoal<Map<String, StreamingResourceStream>, P> fastqMapTransfGoal = new FastqMapTransformGoal(
                getProject(), true, fastqMapGoal, getGoal(GSGoalKey.SETUP));

        FastqDownloadsGoal<P> fastqDownloadsGoal = new FastqDownloadsGoal(getProject(), true, fastqMapGoal, fastqMapTransfGoal,
                getGoal(GSGoalKey.SETUP));

        LoadIndexGoal<P> bloomIndexedGoal = (LoadIndexGoal) getGoal(GSGoalKey.LOAD_INDEX);

        return new FilterGoal<P>(getProject(), fastqMapTransfGoal, bloomIndexedGoal, getExecutionContext(getProject()),
                getGoal(GSGoalKey.SETUP), fastqDownloadsGoal);
    }

    public void filter(boolean clean, String key, String... pathsOrURLs) {
        FilterGoal<P> filterGoal = createGoalChainForFilter(key, pathsOrURLs);
        if (clean) {
            filterGoal.cleanThis();
        }
        filterGoal.make();
    }

    public void filter(String key, String... pathsOrURLs) {
        filter(false, key, pathsOrURLs);
    }

    public void cleanFilter(String key, String... pathsOrURLs) {
        filter(true, key, pathsOrURLs);
    }
}