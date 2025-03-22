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

import org.metagene.genestrip.GSProject.FileType;
import org.metagene.genestrip.bloom.MurmurCGATBloomFilter;
import org.metagene.genestrip.goals.*;
import org.metagene.genestrip.goals.genbank.AssemblyFileDownloadGoal;
import org.metagene.genestrip.goals.genbank.FastaFilesFromGenbankGoal;
import org.metagene.genestrip.goals.genbank.FastaFilesGenbankDownloadGoal;
import org.metagene.genestrip.goals.genbank.TaxNodesFromGenbankGoal;
import org.metagene.genestrip.goals.kraken.KrakenResCountGoal;
import org.metagene.genestrip.goals.refseq.AccessionMapGoal;
import org.metagene.genestrip.goals.refseq.AccessionMapSizeGoal;
import org.metagene.genestrip.goals.refseq.BloomIndexGoal;
import org.metagene.genestrip.goals.refseq.CategoriesGoal;
import org.metagene.genestrip.goals.refseq.CheckSumMapGoal;
import org.metagene.genestrip.goals.refseq.DBGoal;
import org.metagene.genestrip.goals.refseq.FillBloomFilterGoal;
import org.metagene.genestrip.goals.refseq.FillDBGoal;
import org.metagene.genestrip.goals.refseq.FillSizeGoal;
import org.metagene.genestrip.goals.refseq.RefSeqCatalogDownloadGoal;
import org.metagene.genestrip.goals.refseq.RefSeqFnaFilesDownloadGoal;
import org.metagene.genestrip.goals.refseq.RefSeqRNumDownloadGoal;
import org.metagene.genestrip.goals.refseq.StoreDBGoal;
import org.metagene.genestrip.goals.refseq.StoreIndexGoal;
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

public class GSMaker extends Maker<GSProject> {
	private ExecutionContext executionContext;

	public GSMaker(GSProject project) {
		super(project);
	}

	public void dumpAll() {
		dump();
		if (executionContext != null) {
			executionContext.dump();
		}
	}

	protected ExecutionContext createExecutionContext(GSProject project) {
		return new DefaultExecutionContext(project.intConfigValue(GSConfigKey.THREADS),
				project.longConfigValue(GSConfigKey.LOG_PROGRESS_UPDATE_CYCLE));
	}

	protected ExecutionContext getExecutionContext(GSProject project) {
		if (executionContext == null) {
			executionContext = createExecutionContext(project);
		}
		return executionContext;
	}

	protected void createGoals() {
		GSProject project = getProject();
		// Show goals...

		Goal<GSProject> showAllGoals = new Goal<GSProject>(project, GSGoalKey.SHOWALL) {
			@Override
			public boolean isMade() {
				return false;
			}

			@Override
			protected void doMakeThis() {
				List<String> res = new ArrayList<String>();
				for (GSGoalKey key : GSGoalKey.values()) {
					if (!(getGoal(key) instanceof ObjectGoal<?, ?>)) {
						res.add(key.getName());
					}
				}
				System.out.println(res);
			}
		};
		registerGoal(showAllGoals);

		Goal<GSProject> showGoals = new Goal<GSProject>(project, GSGoalKey.SHOW) {
			@Override
			public boolean isMade() {
				return false;
			}

			@Override
			protected void doMakeThis() {
				List<String> res = new ArrayList<String>();
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

		Goal<GSProject> commonSetupGoal = new FileListGoal<GSProject>(project, GSGoalKey.COMMON_SETUP,
				Arrays.asList(project.getCommon().getCommonDir(), project.getCommon().getRefSeqDir(),
						project.getCommon().getGenbankDir(), project.getCommon().getFastqDir(),
						project.getCommon().getFastaDir())) {
			@Override
			protected void makeFile(File file) throws IOException {
				file.mkdir();
			}

			@Override
			protected List<File> getFilesToClean() {
				List<File> files = new ArrayList<File>(getFiles());
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

		Goal<GSProject> taxDBGoal = new TaxIdFileDownloadGoal(project, commonSetupGoal);
		registerGoal(taxDBGoal);

		ObjectGoal<TaxTree, GSProject> taxTreeGoal = new ObjectGoal<TaxTree, GSProject>(project, GSGoalKey.TAXTREE,
				taxDBGoal) {
			@Override
			protected void doMakeThis() {
				set(new TaxTree(getProject().getCommon().getCommonDir()));
			}
		};
		registerGoal(taxTreeGoal);

		ObjectGoal<Set<TaxIdNode>, GSProject> taxNodesGoal = new TaxNodesGoal(project, taxTreeGoal);
		registerGoal(taxNodesGoal);

		// General stuff from RefSeq

		FileGoal<GSProject> releaseNumberGoal = new RefSeqRNumDownloadGoal(project, commonSetupGoal);
		registerGoal(releaseNumberGoal);

		RefSeqCatalogDownloadGoal refSeqCatalogGoal = new RefSeqCatalogDownloadGoal(project, releaseNumberGoal,
				commonSetupGoal);
		registerGoal(refSeqCatalogGoal);

		CheckSumMapGoal checkSumMapGoal = new CheckSumMapGoal(project, refSeqCatalogGoal);
		registerGoal(checkSumMapGoal);

		// Create or clear project directories

		List<File> projectDirs = Arrays.asList(project.getFastaDir(), project.getFastqDir(), project.getDBDir(),
				project.getKrakenOutDir(), project.getResultsDir(), project.getLogDir());

		Goal<GSProject> projectSetupGoal = new FileListGoal<GSProject>(project, GSGoalKey.SETUP, projectDirs,
				commonSetupGoal) {
			@Override
			protected List<File> getFilesToClean() {
				List<File> files = new ArrayList<File>(getFiles());
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

		Goal<GSProject> clearGoal = new FileListGoal<GSProject>(project, GSGoalKey.CLEAR, Arrays
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

		ObjectGoal<Set<RefSeqCategory>, GSProject> categoriesGoal = new CategoriesGoal(project, projectSetupGoal);
		registerGoal(categoriesGoal);

		RefSeqFnaFilesDownloadGoal refSeqFnaFilesGoal = new RefSeqFnaFilesDownloadGoal(project, categoriesGoal,
				refSeqCatalogGoal, checkSumMapGoal);
		registerGoal(refSeqFnaFilesGoal);

		ObjectGoal<Integer, GSProject> accessMapSizeGoal = new AccessionMapSizeGoal(project, categoriesGoal,
				refSeqCatalogGoal);
		registerGoal(accessMapSizeGoal);

		ObjectGoal<AccessionMap, GSProject> accessCollGoal = new AccessionMapGoal(project, categoriesGoal, taxTreeGoal,
				refSeqCatalogGoal, accessMapSizeGoal);
		registerGoal(accessCollGoal);

		// Genbank related additional fastas:
		ObjectGoal<Set<TaxIdNode>, GSProject> taxNodesFromGenBank = new TaxNodesFromGenbankGoal(project, categoriesGoal,
				taxNodesGoal, accessCollGoal);
		registerGoal(taxNodesFromGenBank);

		FileGoal<GSProject> assemblyFileDownloadGoal = new AssemblyFileDownloadGoal(project, commonSetupGoal);
		registerGoal(assemblyFileDownloadGoal);

		FastaFilesFromGenbankGoal fastaFilesFromGenbankGoal = new FastaFilesFromGenbankGoal(project, taxTreeGoal,
				assemblyFileDownloadGoal, taxNodesFromGenBank);
		registerGoal(fastaFilesFromGenbankGoal);

		FastaFilesGenbankDownloadGoal fastaFilesGenbankDownloadGoal = new FastaFilesGenbankDownloadGoal(project,
				fastaFilesFromGenbankGoal);
		registerGoal(fastaFilesGenbankDownloadGoal);
		// end of Genbank stuff.

		AdditionalDownloadsGoal additionalDownloadsGoal = new AdditionalDownloadsGoal(project, projectSetupGoal);
		registerGoal(additionalDownloadsGoal);

		ObjectGoal<Map<File, TaxIdNode>, GSProject> additionalFastasGoal = new AdditionalFastasGoal(project,
				taxTreeGoal, fastaFilesFromGenbankGoal, fastaFilesGenbankDownloadGoal, projectSetupGoal,
				additionalDownloadsGoal);
		registerGoal(additionalFastasGoal);

		// Create database and bloom filter

		FillSizeGoal fillSizeGoal = new FillSizeGoal(project, categoriesGoal, taxNodesGoal, refSeqFnaFilesGoal,
				additionalFastasGoal, accessCollGoal);
		registerGoal(fillSizeGoal);

		ObjectGoal<MurmurCGATBloomFilter, GSProject> fillBloomGoal = new FillBloomFilterGoal(project, categoriesGoal,
				taxNodesGoal, refSeqFnaFilesGoal, additionalFastasGoal, accessCollGoal, fillSizeGoal);
		registerGoal(fillBloomGoal);

		FillDBGoal fillDBGoal = new FillDBGoal(project, categoriesGoal, taxNodesGoal, taxTreeGoal, refSeqFnaFilesGoal,
				additionalFastasGoal, accessCollGoal, fillBloomGoal, projectSetupGoal);
		registerGoal(fillDBGoal);

		StoreDBGoal storeTempDBGoal = new StoreDBGoal(project, GSGoalKey.TEMPDB,
				project.getOutputFile(GSGoalKey.TEMPDB.getName(), FileType.DB, false), fillDBGoal, projectSetupGoal) {
			@Override
			protected void dependentMade(Goal<GSProject> goal) {
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

		FilledDBGoal filledDBGoal = new FilledDBGoal(project, fillDBGoal, storeTempDBGoal);
		registerGoal(filledDBGoal);

		Goal<GSProject> tempDbInfoGoal = new DBInfoGoal(true, project, filledDBGoal, projectSetupGoal);
		registerGoal(tempDbInfoGoal);

		DBGoal updateDBGoal = new DBGoal(project, getExecutionContext(project), categoriesGoal, taxNodesGoal, taxTreeGoal,
				refSeqFnaFilesGoal, additionalFastasGoal, accessCollGoal, filledDBGoal, tempDbInfoGoal, projectSetupGoal);
		registerGoal(updateDBGoal);

		// We weed the storeTempDBGoal as a dependency here so that it does not get
		// automatically deleted until
		// the final database got stored.
		StoreDBGoal storeDBGoal = new StoreDBGoal(project, GSGoalKey.DB, project.getDBFile(), updateDBGoal, projectSetupGoal);
		registerGoal(storeDBGoal);

		LoadDBGoal loadDBGoal = new LoadDBGoal(project, updateDBGoal, storeDBGoal);
		registerGoal(loadDBGoal);

		BloomIndexGoal bloomIndexGoal = new BloomIndexGoal(project, loadDBGoal, projectSetupGoal);
		registerGoal(bloomIndexGoal);

		StoreIndexGoal storeIndexGoal = new StoreIndexGoal(project, bloomIndexGoal, projectSetupGoal);
		registerGoal(storeIndexGoal);

		LoadIndexGoal bloomIndexedGoal = new LoadIndexGoal(project, bloomIndexGoal, storeIndexGoal, projectSetupGoal);
		registerGoal(bloomIndexedGoal);

		// Manage the database

		Goal<GSProject> dbInfoGoal = new DBInfoGoal(project, loadDBGoal, projectSetupGoal);
		registerGoal(dbInfoGoal);

		Goal<GSProject> all = new Goal<GSProject>(project, GSGoalKey.GENALL, dbInfoGoal, storeIndexGoal) {
			@Override
			public boolean isMade() {
				return false;
			}

			@Override
			protected void doMakeThis() {
			}
		};
		registerDefaultGoal(all);

		ObjectGoal<Set<SmallTaxIdNode>, GSProject> db2fastqTaxNodesGoal = new DB2FastqTaxNodesGoal(project, loadDBGoal,
				projectSetupGoal);
		registerGoal(db2fastqTaxNodesGoal);

		Goal<GSProject> db2fastqGoal = new DB2FastqGoal(project, db2fastqTaxNodesGoal, loadDBGoal, projectSetupGoal);
		registerGoal(db2fastqGoal);

		// Use database and bloom filter

		ObjectGoal<Map<String, StreamingResourceStream>, GSProject> fastqMapGoal = new FastqMapGoal(project,
				projectSetupGoal);
		registerGoal(fastqMapGoal);

		ObjectGoal<Map<String, StreamingResourceStream>, GSProject> fastqMapTransfGoal = new FastqMapTransformGoal(
				project, fastqMapGoal, projectSetupGoal);
		registerGoal(fastqMapTransfGoal);

		FastqDownloadsGoal fastqDownloadsGoal = new FastqDownloadsGoal(project, fastqMapGoal, fastqMapTransfGoal,
				projectSetupGoal);
		registerGoal(fastqDownloadsGoal);

		Goal<GSProject> filterGoal = new FilterGoal(project, fastqMapTransfGoal, bloomIndexedGoal,
				getExecutionContext(project), projectSetupGoal, fastqDownloadsGoal);
		registerGoal(filterGoal);

		ObjectGoal<Map<String, MatchingResult>, GSProject> matchResGoal = new MatchResultGoal(getProject(), GSGoalKey.MATCHRES, fastqMapTransfGoal, loadDBGoal,
				getExecutionContext(getProject()), getGoal(GSGoalKey.SETUP), fastqDownloadsGoal);
		registerGoal(matchResGoal);

		Goal<GSProject> matchGoal = new MatchGoal(project, GSGoalKey.MATCH, fastqMapTransfGoal, matchResGoal, projectSetupGoal, fastqDownloadsGoal);
		registerGoal(matchGoal);

		ObjectGoal<Map<String, MatchingResult>, GSProject> matchReslrGoal = new MatchResultGoal(getProject(), GSGoalKey.MATCHRESLR, fastqMapTransfGoal, loadDBGoal,
				getExecutionContext(getProject()), getGoal(GSGoalKey.SETUP), fastqDownloadsGoal);
		registerGoal(matchReslrGoal);

		Goal<GSProject> matchlrGoal = new MatchGoal(project, GSGoalKey.MATCHLR, fastqMapTransfGoal, matchReslrGoal, projectSetupGoal, fastqDownloadsGoal);
		registerGoal(matchlrGoal);

		// Use kraken
		KrakenResCountGoal krakenResCountGoal = new KrakenResCountGoal(project, fastqMapTransfGoal, taxNodesGoal,
				projectSetupGoal, fastqDownloadsGoal);
		registerGoal(krakenResCountGoal);
	}

	public MatchingResult cleanMatch(boolean lr, String key, String... pathsOrURLs) {
		return match(lr, true, key, pathsOrURLs);
	}

	public MatchingResult match(boolean lr, String key, String... pathsOrURLs) {
		return match(lr, false, key, pathsOrURLs);
	}

	public MatchingResult match(boolean lr, boolean clean, String key, String... pathsOrURLs) {
		MatchResultGoal matchResGoal = createGoalChainForMatchResult(lr, key, pathsOrURLs);
		LoadDBGoal loadDBGoal = (LoadDBGoal) getGoal(GSGoalKey.LOAD_DB);
		MatchGoal matchGoal = new MatchGoal(matchResGoal.getProject(), (lr ? GSGoalKey.MATCHLR : GSGoalKey.MATCH), matchResGoal.getFastqMapGoal(), matchResGoal);
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

	protected MatchResultGoal createGoalChainForMatchResult(boolean lr, String key, String... pathsOrURLs) {
		ObjectGoal<Map<String, StreamingResourceStream>, GSProject> fastqMapGoal = new FastqMapGoal(getProject(),
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

		ObjectGoal<Map<String, StreamingResourceStream>, GSProject> fastqMapTransfGoal = new FastqMapTransformGoal(
				getProject(), fastqMapGoal, getGoal(GSGoalKey.SETUP));

		FastqDownloadsGoal fastqDownloadsGoal = new FastqDownloadsGoal(getProject(), fastqMapGoal, fastqMapTransfGoal,
				getGoal(GSGoalKey.SETUP));

		LoadDBGoal loadDBGoal = (LoadDBGoal) getGoal(GSGoalKey.LOAD_DB);

		return new MatchResultGoal(getProject(), (lr ? GSGoalKey.MATCHRESLR : GSGoalKey.MATCHRES), fastqMapTransfGoal, loadDBGoal,
				getExecutionContext(getProject()), getGoal(GSGoalKey.SETUP), fastqDownloadsGoal);
/*		return new MatchResultGoal(getProject(), (lr ? GSGoalKey.MATCHRESLR : GSGoalKey.MATCHRES), fastqMapTransfGoal, loadDBGoal,
				getExecutionContext(getProject()), getGoal(GSGoalKey.SETUP), fastqDownloadsGoal);
/*		{
			protected void writeOutputFile(File file, MatchingResult result, Database wrapper) throws IOException {
				if (key == null || !key.isEmpty()) {
					super.writeOutputFile(file, result);
				}
			}
		};*/
	}

	protected FilterGoal createGoalChainForFilter(String key, String... pathsOrURLs) {
		ObjectGoal<Map<String, StreamingResourceStream>, GSProject> fastqMapGoal = new FastqMapGoal(getProject(),
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

		ObjectGoal<Map<String, StreamingResourceStream>, GSProject> fastqMapTransfGoal = new FastqMapTransformGoal(
				getProject(), fastqMapGoal, getGoal(GSGoalKey.SETUP));

		FastqDownloadsGoal fastqDownloadsGoal = new FastqDownloadsGoal(getProject(), fastqMapGoal, fastqMapTransfGoal,
				getGoal(GSGoalKey.SETUP));

		LoadIndexGoal bloomIndexedGoal = (LoadIndexGoal) getGoal(GSGoalKey.LOAD_INDEX);

		return new FilterGoal(getProject(), fastqMapTransfGoal, bloomIndexedGoal, getExecutionContext(getProject()),
				getGoal(GSGoalKey.SETUP), fastqDownloadsGoal);
	}

	public void filter(boolean clean, String key, String... pathsOrURLs) {
		FilterGoal filterGoal = createGoalChainForFilter(key, pathsOrURLs);
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