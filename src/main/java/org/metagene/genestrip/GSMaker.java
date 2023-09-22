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
import org.metagene.genestrip.goals.AssemblyFileDownloadGoal;
import org.metagene.genestrip.goals.BloomFilterFileGoal;
import org.metagene.genestrip.goals.BloomFilterSizeGoal;
import org.metagene.genestrip.goals.FastaFileDownloadGoal;
import org.metagene.genestrip.goals.FilterGoal;
import org.metagene.genestrip.goals.KMerFastqGoal;
import org.metagene.genestrip.goals.KMerFastqStoreFileGoal;
import org.metagene.genestrip.goals.KMerStoreFileGoal;
import org.metagene.genestrip.goals.KrakenFastqFileGoal;
import org.metagene.genestrip.goals.KrakenOutGoal;
import org.metagene.genestrip.goals.KrakenResCountGoal;
import org.metagene.genestrip.goals.KrakenResErrorGoal;
import org.metagene.genestrip.goals.MatchGoal;
import org.metagene.genestrip.goals.MultiMatchGoal;
import org.metagene.genestrip.goals.SortKrakenOutGoal;
import org.metagene.genestrip.goals.StoreInfoGoal;
import org.metagene.genestrip.goals.TaxIdFileDownloadGoal;
import org.metagene.genestrip.goals.TaxNodesGoal;
import org.metagene.genestrip.goals.TrieFromKrakenResGoal;
import org.metagene.genestrip.goals.TrieFromKrakenResGoal.TaxidWithCount;
import org.metagene.genestrip.goals.refseq.AccessionMap;
import org.metagene.genestrip.goals.refseq.AccessionMapGoal;
import org.metagene.genestrip.goals.refseq.AccessionMapSizeGoal;
import org.metagene.genestrip.goals.refseq.CategoriesGoal;
import org.metagene.genestrip.goals.refseq.FillBloomFilterGoal;
import org.metagene.genestrip.goals.refseq.FillSizeGoal;
import org.metagene.genestrip.goals.refseq.FillStoreGoal;
import org.metagene.genestrip.goals.refseq.RefSeqCatalogDownloadGoal;
import org.metagene.genestrip.goals.refseq.RefSeqCategory;
import org.metagene.genestrip.goals.refseq.RefSeqFnaFilesDownloadGoal;
import org.metagene.genestrip.goals.refseq.RefSeqRNumDownloadGoal;
import org.metagene.genestrip.goals.refseq.UpdateStoreGoal;
import org.metagene.genestrip.make.FileGoal;
import org.metagene.genestrip.make.FileListGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.Maker;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.store.KMerTrie;
import org.metagene.genestrip.tax.AdditionalFastaInfoReader;
import org.metagene.genestrip.tax.AssemblySummaryReader;
import org.metagene.genestrip.tax.AssemblySummaryReader.FTPEntryWithQuality;
import org.metagene.genestrip.tax.TaxTree;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;

public class GSMaker extends Maker<GSProject> {
	public GSMaker(GSProject project) {
		super(project);
	}

	protected void createGoals(GSProject project) {
		List<File> projectDirs = Arrays.asList(project.getFastasDir(), project.getFastqsDir(), project.getFiltersDir(),
				project.getKrakenOutDir(), project.getResultsDir());

		Goal<GSProject> commonSetupGoal = new FileListGoal<GSProject>(project, "commonsetup",
				Arrays.asList(project.getConfig().getCommonDir(), project.getConfig().getRefSeqDir())) {
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

		Goal<GSProject> clearGoal = new FileListGoal<GSProject>(project, "clear", Arrays.asList(project.getFastqsDir(),
				project.getFiltersDir(), project.getKrakenOutDir(), project.getResultsDir())) {
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

		FileGoal<GSProject> assemblyGoal = new AssemblyFileDownloadGoal(project, "assemblydownload", commonSetupGoal);
		registerGoal(assemblyGoal);

		ObjectGoal<Map<TaxIdNode, List<FTPEntryWithQuality>>, GSProject> fastaFilesGoal = new ObjectGoal<Map<TaxIdNode, List<FTPEntryWithQuality>>, GSProject>(
				project, "fastafiles", assemblyGoal, taxNodesGoal) {
			@Override
			public void makeThis() {
				try {
					AssemblySummaryReader assemblySummaryReader = new AssemblySummaryReader(
							getProject().getConfig().getCommonDir(), getProject().getConfig().isUseGenBank(),
							taxTreeGoal.get());
					int[] nEntriesTotal = new int[1];
					Map<TaxIdNode, List<FTPEntryWithQuality>> entries = assemblySummaryReader.getRelevantEntries(
							taxNodesGoal.get(), getProject().getConfig().getFastaQualities(), nEntriesTotal);
					if (getLogger().isInfoEnabled()) {
						getLogger().info("Total number of entries in assembly summary file: " + nEntriesTotal[0]);
					}

					AdditionalFastaInfoReader additionalFastaInfoReader = new AdditionalFastaInfoReader(
							getProject().getConfig().getAdditionalDir(), taxTreeGoal.get());
					additionalFastaInfoReader.addRelevantEntries(entries, taxNodesGoal.get(), nEntriesTotal);
					if (getLogger().isInfoEnabled()) {
						getLogger().info("Total number of entries in additonal info file: " + nEntriesTotal[0]);
					}

					set(entries);
				} catch (IOException e) {
					throw new RuntimeException(e);
				}
			}
		};

		FastaFileDownloadGoal fastaDownloadGoal = new FastaFileDownloadGoal(project, "fastasdownload", fastaFilesGoal,
				projectSetupGoal);
		registerGoal(fastaDownloadGoal);

		KMerFastqGoal kmerFastqGoal = new KMerFastqGoal(project, "kmerfastq", fastaFilesGoal, fastaDownloadGoal,
				projectSetupGoal);
		registerGoal(kmerFastqGoal);

		FileGoal<GSProject> krakenOutGoal = new KrakenOutGoal(project, "kmerkrakenout", kmerFastqGoal,
				projectSetupGoal);
		registerGoal(krakenOutGoal);

		if (project.getConfig().isUseKraken1()) {
			krakenOutGoal = new SortKrakenOutGoal(project, "sort", krakenOutGoal);
			registerGoal(krakenOutGoal);
		}

		BloomFilterSizeGoal bloomFilterSizeGoal = new BloomFilterSizeGoal(project, "bloomsize", taxNodesGoal,
				krakenOutGoal);
		registerGoal(bloomFilterSizeGoal);

		KMerStoreFileGoal storeGoal = new KMerStoreFileGoal(project, "oldstore", taxTreeGoal, taxNodesGoal,
				krakenOutGoal, kmerFastqGoal, bloomFilterSizeGoal, projectSetupGoal, taxTreeGoal);
		registerGoal(storeGoal);

		KrakenFastqFileGoal krakenFastqGoal = new KrakenFastqFileGoal(project, "krakenfastq", taxNodesGoal,
				krakenOutGoal, kmerFastqGoal, projectSetupGoal);
		registerGoal(krakenFastqGoal);

		Goal<GSProject> kMerFastqStoreFileGoal = new KMerFastqStoreFileGoal(project, "store2", taxNodesGoal,
				krakenFastqGoal, projectSetupGoal);
		registerGoal(kMerFastqStoreFileGoal);

		BloomFilterFileGoal bloomFilterFileGoal = new BloomFilterFileGoal(project, "bloom", bloomFilterSizeGoal,
				taxNodesGoal, krakenOutGoal, kmerFastqGoal, projectSetupGoal);
		registerGoal(bloomFilterFileGoal);

		Goal<GSProject> showGoals = new Goal<GSProject>(project, "show") {
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
		registerGoal(showGoals);

		Goal<GSProject> all = new Goal<GSProject>(project, "genall", kMerFastqStoreFileGoal, bloomFilterFileGoal) {
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
			Goal<GSProject> filterGoal = new FilterGoal(project, "filter", fastqOrCSV,
					project.getConfig().isWriteDumpedFastq(), bloomFilterFileGoal, projectSetupGoal);
			registerGoal(filterGoal);

			Goal<GSProject> classifyGoal = new MatchGoal(project, "oldmatch", fastqOrCSV, taxTreeGoal, storeGoal,
					project.getConfig().isWriteFilteredFastq(), projectSetupGoal);
			registerGoal(classifyGoal);

			Goal<GSProject> multiMatchGoal = new MultiMatchGoal(project, "old" + MultiMatchGoal.NAME, fastqOrCSV, taxTreeGoal,
					storeGoal, project.getConfig().isWriteFilteredFastq(), projectSetupGoal);
			registerGoal(multiMatchGoal);

			Goal<GSProject> krakenResCountGoal = new KrakenResCountGoal(project, "krakenres", fastqOrCSV, taxNodesGoal,
					projectSetupGoal);
			registerGoal(krakenResCountGoal);

			Goal<GSProject> multiKrakenResCountGoal = new KrakenResCountGoal(project, "multikrakenres", fastqOrCSV,
					true, taxNodesGoal, projectSetupGoal);
			registerGoal(multiKrakenResCountGoal);

			Goal<GSProject> krakenResCountAllGoal = new KrakenResCountGoal(project, "krakenresall", fastqOrCSV, null,
					projectSetupGoal);
			registerGoal(krakenResCountAllGoal);

			Goal<GSProject> multiKrakenResCountAllGoal = new KrakenResCountGoal(project, "multikrakenresall",
					fastqOrCSV, true, null, projectSetupGoal);
			registerGoal(multiKrakenResCountAllGoal);

			FileGoal<GSProject> fastqKrakenOutGoal = new KrakenOutGoal(project, "fastqkrakenout", fastqOrCSV,
					projectSetupGoal);
			registerGoal(fastqKrakenOutGoal);

			if (project.getConfig().isUseKraken1()) {
				fastqKrakenOutGoal = new SortKrakenOutGoal(project, "sortfastq", fastqKrakenOutGoal);
				registerGoal(fastqKrakenOutGoal);
			}

			ObjectGoal<KMerTrie<TaxidWithCount>, GSProject> trieFromKrakenResGoal = new TrieFromKrakenResGoal(project,
					"triefromkrakenres", fastqOrCSV, project.getConfig().isWriteFilteredFastq(), taxNodesGoal,
					fastaFilesGoal, fastqKrakenOutGoal, fastaDownloadGoal, projectSetupGoal);
			registerGoal(trieFromKrakenResGoal);

			Goal<GSProject> krakenResErrorGoal = new KrakenResErrorGoal(project, "krakenerr", fastqOrCSV,
					trieFromKrakenResGoal, projectSetupGoal);
			registerGoal(krakenResErrorGoal);
		}

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

		FillSizeGoal fillSizeGoal = new FillSizeGoal(project, "fillsize", categoriesGoal, taxNodesGoal,
				refSeqFnaFilesGoal, accessCollGoal);
		registerGoal(fillSizeGoal);

		ObjectGoal<MurmurCGATBloomFilter, GSProject> fillBloomGoal = new FillBloomFilterGoal(project, "fillbloom",
				categoriesGoal, taxNodesGoal, refSeqFnaFilesGoal, accessCollGoal, fillSizeGoal);
		registerGoal(fillBloomGoal);

		FillStoreGoal fillStoreGoal = new FillStoreGoal(project, "tempstore", categoriesGoal, taxNodesGoal,
				refSeqFnaFilesGoal, accessCollGoal, fillSizeGoal, fillBloomGoal, projectSetupGoal);
		registerGoal(fillStoreGoal);

		FileGoal<GSProject> updateStoreGoal = new UpdateStoreGoal(project, "store", categoriesGoal, taxTreeGoal,
				refSeqFnaFilesGoal, accessCollGoal, fillStoreGoal, projectSetupGoal);
		registerGoal(updateStoreGoal);

		Goal<GSProject> newMatchGoal = new MatchGoal(project, "match", fastqOrCSV, taxTreeGoal, updateStoreGoal,
				project.getConfig().isWriteFilteredFastq(), projectSetupGoal);
		registerGoal(newMatchGoal);

		Goal<GSProject> newMultiMatchGoal = new MultiMatchGoal(project, MultiMatchGoal.NAME, fastqOrCSV, taxTreeGoal,
				updateStoreGoal, project.getConfig().isWriteFilteredFastq(), projectSetupGoal);
		registerGoal(newMultiMatchGoal);

		Goal<GSProject> storeInfoGoal = new StoreInfoGoal(project, "storeinfo", taxTreeGoal, updateStoreGoal,
				projectSetupGoal);
		registerGoal(storeInfoGoal);
	}
}