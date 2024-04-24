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

import org.metagene.genestrip.bloom.MurmurCGATBloomFilter;
import org.metagene.genestrip.goals.AdditionalFastasGoal;
import org.metagene.genestrip.goals.DB2FastqGoal;
import org.metagene.genestrip.goals.DB2FastqTaxNodesGoal;
import org.metagene.genestrip.goals.FilterGoal;
import org.metagene.genestrip.goals.MultiMatchGoal;
import org.metagene.genestrip.goals.StoreInfoGoal;
import org.metagene.genestrip.goals.TaxIdFileDownloadGoal;
import org.metagene.genestrip.goals.TaxNodesGoal;
import org.metagene.genestrip.goals.genbank.AssemblyFileDownloadGoal;
import org.metagene.genestrip.goals.genbank.FastaFilesFromGenbankGoal;
import org.metagene.genestrip.goals.genbank.FastaFilesGenbankDownloadGoal;
import org.metagene.genestrip.goals.genbank.TaxNodesFromGenbankGoal;
import org.metagene.genestrip.goals.kraken.KrakenResCountGoal;
import org.metagene.genestrip.goals.refseq.AccessionMapGoal;
import org.metagene.genestrip.goals.refseq.AccessionMapSizeGoal;
import org.metagene.genestrip.goals.refseq.BloomIndexGoal;
import org.metagene.genestrip.goals.refseq.BloomIndexedGoal;
import org.metagene.genestrip.goals.refseq.CategoriesGoal;
import org.metagene.genestrip.goals.refseq.FillBloomFilterGoal;
import org.metagene.genestrip.goals.refseq.FillSizeGoal;
import org.metagene.genestrip.goals.refseq.FillStoreGoal;
import org.metagene.genestrip.goals.refseq.FilledStoreGoal;
import org.metagene.genestrip.goals.refseq.RefSeqCatalogDownloadGoal;
import org.metagene.genestrip.goals.refseq.RefSeqFnaFilesDownloadGoal;
import org.metagene.genestrip.goals.refseq.RefSeqRNumDownloadGoal;
import org.metagene.genestrip.goals.refseq.UpdateStoreGoal;
import org.metagene.genestrip.goals.refseq.UpdatedStoreGoal;
import org.metagene.genestrip.make.FileGoal;
import org.metagene.genestrip.make.FileListGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.Maker;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.refseq.AccessionMap;
import org.metagene.genestrip.refseq.RefSeqCategory;
import org.metagene.genestrip.tax.TaxTree;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;
import org.metagene.genestrip.util.DefaultExecutorServiceBundle;
import org.metagene.genestrip.util.ExecutorServiceBundle;

public class GSMaker extends Maker<GSProject> {
	public enum UserGoal {
		SHOW("show"), SHOWALL("showall"), DB("db"), DBINFO("dbinfo"), DB2FASTQ("db2fastq"), INDEX("index"),
		GENALL("genall"), CLEAR("clear"), MATCH("match"), MULTIMATCH("multimatch"), FILTER("filtler"),
		MATCHLR("matchlr"), MULTIMATCHLR("multimatchlr");

		private final String name;

		private UserGoal(String name) {
			this.name = name;
		}

		public String getName() {
			return name;
		}
	}

	private ExecutorServiceBundle executorServiceBundle;

	public GSMaker(GSProject project) {
		super(project);
	}

	public Goal<GSProject> getGoal(UserGoal goal) {
		return getGoal(goal.getName());
	}

	public void dump() {
		if (executorServiceBundle != null) {
			executorServiceBundle.dump();
		}
	}

	protected ExecutorServiceBundle createExecutorServiceBundle(GSConfig config) {
		return new DefaultExecutorServiceBundle(config.getThreads(), config.getMatchLogUpdateCycle()) {
			protected String getThreadBaseName() {
				return "Fastq Reader Thread";
			}
		};
	}

	protected ExecutorServiceBundle getExecutorServiceBundle(GSConfig config) {
		if (executorServiceBundle == null) {
			executorServiceBundle = createExecutorServiceBundle(config);
		}
		return executorServiceBundle;
	}

	protected void createGoals(GSProject project) {
		List<File> projectDirs = Arrays.asList(project.getFastaDir(), project.getFastqDir(), project.getDBDir(),
				project.getKrakenOutDir(), project.getResultsDir(), project.getLogDir(), project.getGenbankDir());

		Goal<GSProject> commonSetupGoal = new FileListGoal<GSProject>(project, "commonsetup",
				Arrays.asList(project.getConfig().getCommonDir(), project.getConfig().getRefSeqDir(),
						project.getConfig().getGenbankDir(), project.getConfig().getFastqDir(),
						project.getConfig().getFastaDir())) {
			@Override
			protected void makeFile(File file) throws IOException {
				file.mkdir();
			}

			@Override
			public boolean isAllowTransitiveClean() {
				return false;
			}
		};
		registerGoal(commonSetupGoal);

		Goal<GSProject> projectSetupGoal = new FileListGoal<GSProject>(project, "setup", projectDirs) {
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

		Goal<GSProject> clearGoal = new FileListGoal<GSProject>(project, UserGoal.CLEAR.getName(),
				Arrays.asList(project.getDBDir(), project.getKrakenOutDir(), project.getResultsDir())) {
			@Override
			public boolean isMade() {
				return false;
			}

			@Override
			protected void makeFile(File file) throws IOException {
			}

			@Override
			public void makeThis() {
				cleanThis();
			}
		};
		registerGoal(clearGoal);

		Goal<GSProject> taxDBGoal = new TaxIdFileDownloadGoal(project, "taxdownload", commonSetupGoal);
		registerGoal(taxDBGoal);

		ObjectGoal<TaxTree, GSProject> taxTreeGoal = new ObjectGoal<TaxTree, GSProject>(project, "taxtree", taxDBGoal) {
			@Override
			public void makeThis() {
				set(new TaxTree(getProject().getConfig().getCommonDir()));
			}
		};
		registerGoal(taxTreeGoal);

		ObjectGoal<Set<TaxIdNode>, GSProject> taxNodesGoal = new TaxNodesGoal(project, "taxids", taxTreeGoal);
		registerGoal(taxNodesGoal);

		Goal<GSProject> showAllGoals = new Goal<GSProject>(project, UserGoal.SHOWALL.getName()) {
			@Override
			public boolean isMade() {
				return false;
			}

			@Override
			public void makeThis() {
				List<String> res = new ArrayList<String>();
				for (String name : getGoalNames()) {
					if (!(getGoal(name) instanceof ObjectGoal<?, ?>)) {
						res.add(name);
					}
				}
				System.out.println(res);
			}
		};
		registerGoal(showAllGoals);

		Goal<GSProject> showGoals = new Goal<GSProject>(project, UserGoal.SHOW.getName()) {
			@Override
			public boolean isMade() {
				return false;
			}

			@Override
			public void makeThis() {
				List<String> res = new ArrayList<String>();
				for (UserGoal goal : UserGoal.values()) {
					res.add(goal.getName());
				}
				System.out.println(res);
			}
		};
		registerGoal(showGoals);

		ObjectGoal<Set<RefSeqCategory>[], GSProject> categoriesGoal = new CategoriesGoal(project, "refseqcats",
				projectSetupGoal);

		FileGoal<GSProject> releaseNumberGoal = new RefSeqRNumDownloadGoal(project, "refseqrelease", commonSetupGoal);
		registerGoal(releaseNumberGoal);

		RefSeqCatalogDownloadGoal refSeqCatalogGoal = new RefSeqCatalogDownloadGoal(project, "refseqcat",
				releaseNumberGoal, commonSetupGoal);
		registerGoal(refSeqCatalogGoal);

		RefSeqFnaFilesDownloadGoal refSeqFnaFilesGoal = new RefSeqFnaFilesDownloadGoal(project, "refseqfna",
				categoriesGoal, refSeqCatalogGoal);
		registerGoal(refSeqFnaFilesGoal);

		ObjectGoal<Integer, GSProject> accessMapSizeGoal = new AccessionMapSizeGoal(project, "accmapsize",
				categoriesGoal, refSeqCatalogGoal, refSeqFnaFilesGoal);
		registerGoal(accessMapSizeGoal);

		ObjectGoal<AccessionMap, GSProject> accessCollGoal = new AccessionMapGoal(project, "accmap", categoriesGoal,
				taxTreeGoal, refSeqCatalogGoal, refSeqFnaFilesGoal, accessMapSizeGoal);
		registerGoal(accessCollGoal);

		// Genbank related additional fastas:
		ObjectGoal<Set<TaxIdNode>, GSProject> taxNodesFromGenBank = new TaxNodesFromGenbankGoal(project,
				"taxfromgenbank", categoriesGoal, taxNodesGoal, refSeqFnaFilesGoal, accessCollGoal);
		registerGoal(taxNodesFromGenBank);

		FileGoal<GSProject> assemblyFileDownloadGoal = new AssemblyFileDownloadGoal(project, "assemblydownload",
				commonSetupGoal);
		registerGoal(assemblyFileDownloadGoal);

		FastaFilesFromGenbankGoal fastaFilesFromGenbankGoal = new FastaFilesFromGenbankGoal(project, "fastasgenbank",
				taxTreeGoal, assemblyFileDownloadGoal, taxNodesFromGenBank);
		registerGoal(fastaFilesFromGenbankGoal);

		FastaFilesGenbankDownloadGoal fastaFilesGenbankDownloadGoal = new FastaFilesGenbankDownloadGoal(project,
				"fastasgenbankdl", fastaFilesFromGenbankGoal);
		registerGoal(fastaFilesGenbankDownloadGoal);
		// end of Genbank stuff.

		ObjectGoal<Map<File, TaxIdNode>, GSProject> additionalFastasGoal = new AdditionalFastasGoal(project,
				"addfastas", taxTreeGoal, fastaFilesFromGenbankGoal, fastaFilesGenbankDownloadGoal, commonSetupGoal);
		registerGoal(additionalFastasGoal);

		FillSizeGoal fillSizeGoal = new FillSizeGoal(project, "fillsize", categoriesGoal, taxNodesGoal,
				refSeqFnaFilesGoal, additionalFastasGoal, accessCollGoal);
		registerGoal(fillSizeGoal);

		ObjectGoal<MurmurCGATBloomFilter, GSProject> fillBloomGoal = new FillBloomFilterGoal(project, "fillindex",
				categoriesGoal, taxNodesGoal, refSeqFnaFilesGoal, additionalFastasGoal, accessCollGoal, fillSizeGoal);
		registerGoal(fillBloomGoal);

		FillStoreGoal fillStoreGoal = new FillStoreGoal(project, "tempdb", categoriesGoal, taxNodesGoal,
				refSeqFnaFilesGoal, additionalFastasGoal, accessCollGoal, fillSizeGoal, fillBloomGoal,
				projectSetupGoal);
		registerGoal(fillStoreGoal);

		FilledStoreGoal filledStoreGoal = new FilledStoreGoal(project, "filleddb", fillStoreGoal);
		registerGoal(filledStoreGoal);
		fillStoreGoal.setFilledStoreGoal(filledStoreGoal);

		UpdateStoreGoal updateStoreGoal = new UpdateStoreGoal(project, UserGoal.DB.getName(),
				getExecutorServiceBundle(project.getConfig()), categoriesGoal, taxTreeGoal, refSeqFnaFilesGoal,
				additionalFastasGoal, accessCollGoal, filledStoreGoal, projectSetupGoal);
		registerGoal(updateStoreGoal);

		UpdatedStoreGoal updatedStoreGoal = new UpdatedStoreGoal(project, "updateddb", updateStoreGoal);
		registerGoal(updatedStoreGoal);
		updateStoreGoal.setUpdatedStoreGoal(updatedStoreGoal);

		BloomIndexGoal bloomIndexGoal = new BloomIndexGoal(project, UserGoal.INDEX.getName(), taxTreeGoal, taxNodesGoal,
				updatedStoreGoal, projectSetupGoal);
		registerGoal(bloomIndexGoal);

		BloomIndexedGoal bloomIndexedGoal = new BloomIndexedGoal(project, "indexed", bloomIndexGoal, projectSetupGoal);
		registerGoal(bloomIndexedGoal);
		bloomIndexGoal.setBloomIndexedGoal(bloomIndexedGoal);

		Goal<GSProject> storeInfoGoal = new StoreInfoGoal(project, UserGoal.DBINFO.getName(), taxTreeGoal,
				updatedStoreGoal, projectSetupGoal);
		registerGoal(storeInfoGoal);

		ObjectGoal<Set<TaxIdNode>, GSProject> db2fastqTaxNodesGoal = new DB2FastqTaxNodesGoal(project, "db2fastqtaxids",
				taxTreeGoal, updatedStoreGoal, projectSetupGoal);
		registerGoal(db2fastqTaxNodesGoal);

		Goal<GSProject> db2fastqGoal = new DB2FastqGoal(project, UserGoal.DB2FASTQ.getName(), db2fastqTaxNodesGoal,
				updatedStoreGoal, projectSetupGoal);
		registerGoal(db2fastqGoal);

		Goal<GSProject> all = new Goal<GSProject>(project, UserGoal.GENALL.getName(), storeInfoGoal, bloomIndexGoal) {
			@Override
			public boolean isMade() {
				return false;
			}

			@Override
			public void makeThis() {
			}
		};
		registerDefaultGoal(all);

		File fastqOrCSV = project.getFastqOrCSVFile();
		if (fastqOrCSV != null) {
			Goal<GSProject> filterGoal = new FilterGoal(project, UserGoal.FILTER.getName(), fastqOrCSV,
					bloomIndexedGoal, getExecutorServiceBundle(project.getConfig()), projectSetupGoal);
			registerGoal(filterGoal);

			KrakenResCountGoal krakenResCountGoal = new KrakenResCountGoal(project, "krakenres", false, fastqOrCSV,
					taxNodesGoal, projectSetupGoal);
			registerGoal(krakenResCountGoal);

			KrakenResCountGoal multiKrakenResCountGoal = new KrakenResCountGoal(project, "multikrakenres", true,
					fastqOrCSV, taxNodesGoal, projectSetupGoal);
			registerGoal(multiKrakenResCountGoal);

			Goal<GSProject> krakenResCountAllGoal = new KrakenResCountGoal(project, "krakenresall", false, fastqOrCSV,
					null, projectSetupGoal);
			registerGoal(krakenResCountAllGoal);

			Goal<GSProject> multiKrakenResCountAllGoal = new KrakenResCountGoal(project, "multikrakenresall", true,
					fastqOrCSV, null, projectSetupGoal);
			registerGoal(multiKrakenResCountAllGoal);

			Goal<GSProject> matchGoal = new MultiMatchGoal(project, UserGoal.MATCH.getName(), false, fastqOrCSV,
					taxTreeGoal, updatedStoreGoal, project.getConfig().isWriteFilteredFastq(),
					getExecutorServiceBundle(project.getConfig()), projectSetupGoal);
			registerGoal(matchGoal);

			Goal<GSProject> matchlrGoal = new MultiMatchGoal(project, UserGoal.MATCHLR.getName(), false, fastqOrCSV,
					taxTreeGoal, updatedStoreGoal, project.getConfig().isWriteFilteredFastq(),
					getExecutorServiceBundle(project.getConfig()), projectSetupGoal);
			registerGoal(matchlrGoal);
		}
		if (fastqOrCSV != null || project.getKeyToFastqs() != null) {
			Goal<GSProject> multiMatchGoal = new MultiMatchGoal(project, UserGoal.MULTIMATCH.getName(), true,
					fastqOrCSV, taxTreeGoal, updatedStoreGoal, project.getConfig().isWriteFilteredFastq(),
					getExecutorServiceBundle(project.getConfig()), projectSetupGoal);
			registerGoal(multiMatchGoal);

			Goal<GSProject> multiMatchlrGoal = new MultiMatchGoal(project, UserGoal.MULTIMATCHLR.getName(), true,
					fastqOrCSV, taxTreeGoal, updatedStoreGoal, project.getConfig().isWriteFilteredFastq(),
					getExecutorServiceBundle(project.getConfig()), projectSetupGoal);
			registerGoal(multiMatchlrGoal);
		}
	}
}