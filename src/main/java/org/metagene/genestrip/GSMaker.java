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
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.metagene.genestrip.goals.AssemblyFileDownloadGoal;
import org.metagene.genestrip.goals.BloomFilterFileGoal;
import org.metagene.genestrip.goals.BloomFilterSizeGoal;
import org.metagene.genestrip.goals.ClassifyGoal;
import org.metagene.genestrip.goals.FastaFileDownloadGoal;
import org.metagene.genestrip.goals.FastasSizeGoal;
import org.metagene.genestrip.goals.FilterGoal;
import org.metagene.genestrip.goals.KMerFastqGoal;
import org.metagene.genestrip.goals.KMerFastqTrieFileGoal;
import org.metagene.genestrip.goals.KMerTrieFileGoal;
import org.metagene.genestrip.goals.KrakenFastqFileGoal;
import org.metagene.genestrip.goals.KrakenOutGoal;
import org.metagene.genestrip.goals.KrakenResCountGoal;
import org.metagene.genestrip.goals.KrakenResErrorGoal;
import org.metagene.genestrip.goals.TaxIdFileDownloadGoal;
import org.metagene.genestrip.goals.TaxNodesGoal;
import org.metagene.genestrip.goals.TrieFromKrakenResGoal;
import org.metagene.genestrip.goals.TrieFromKrakenResGoal.TaxidWithCount;
import org.metagene.genestrip.make.FileGoal;
import org.metagene.genestrip.make.FileListGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.Maker;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.tax.AssemblySummaryReader;
import org.metagene.genestrip.tax.AssemblySummaryReader.FTPEntryWithQuality;
import org.metagene.genestrip.tax.TaxTree;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;
import org.metagene.genestrip.trie.KMerTrie;

public class GSMaker extends Maker<GSProject> {
	public GSMaker(GSProject project) {
		super(project);
	}

	protected void createGoals(GSProject project) {
		List<File> projectDirs = Arrays.asList(project.getFastasDir(), project.getFastqsDir(), project.getFiltersDir(),
				project.getKrakenOutDir(), project.getResultsDir());

		Goal<GSProject> commonSetupGoal = new FileListGoal<GSProject>(project, "commonsetup",
				project.getConfig().getCommonDir()) {
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

		Goal<GSProject> commonGoal = new Goal<GSProject>(project, "common", assemblyGoal, taxDBGoal) {
			@Override
			public boolean isMade() {
				return false;
			}
		};
		registerGoal(commonGoal);

		ObjectGoal<Map<TaxIdNode, List<FTPEntryWithQuality>>, GSProject> fastaFilesGoal = new ObjectGoal<Map<TaxIdNode, List<FTPEntryWithQuality>>, GSProject>(
				project, "fastafiles", assemblyGoal, taxNodesGoal) {
			@Override
			public void makeThis() {
				try {
					AssemblySummaryReader assemblySummaryReader = new AssemblySummaryReader(
							getProject().getConfig().getCommonDir(), getProject().getConfig().isUseGenBank(),
							taxTreeGoal.get());
					int[] nEntriesTotal = new int[1];
					Map<TaxIdNode, List<FTPEntryWithQuality>> entries = assemblySummaryReader
							.getRelevantEntries(taxNodesGoal.get(), nEntriesTotal);
					if (getLogger().isInfoEnabled()) {
						getLogger().info("Total number of entries in assembly summary file: " + nEntriesTotal[0]);
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

		FastasSizeGoal fastasSizeGoal = new FastasSizeGoal(project, "fastassize", 10, fastaDownloadGoal);

		KMerFastqGoal kmerFastqGoal = new KMerFastqGoal(project, "kmerfastqgen", fastasSizeGoal, fastaDownloadGoal,
				projectSetupGoal);
		registerGoal(kmerFastqGoal);

		KrakenOutGoal krakenOutGoal = new KrakenOutGoal(project, "kmerkrakenout", kmerFastqGoal.getOutputFile(),
				kmerFastqGoal, projectSetupGoal);
		registerGoal(krakenOutGoal);

		KMerTrieFileGoal trieGoal = new KMerTrieFileGoal(project, "trie", taxNodesGoal, krakenOutGoal, kmerFastqGoal,
				projectSetupGoal);
		registerGoal(trieGoal);

		KrakenFastqFileGoal krakenFastqGoal = new KrakenFastqFileGoal(project, "krakenfastq", taxNodesGoal,
				krakenOutGoal, kmerFastqGoal, projectSetupGoal);
		registerGoal(krakenFastqGoal);

		Goal<GSProject> kMerFastqTrieFileGoal = new KMerFastqTrieFileGoal(project, "trieg2", taxNodesGoal,
				krakenFastqGoal, projectSetupGoal);
		registerGoal(kMerFastqTrieFileGoal);

		BloomFilterSizeGoal bloomFilterSizeGoal = new BloomFilterSizeGoal(project, "bloomsize", taxNodesGoal,
				krakenOutGoal);
		registerGoal(bloomFilterSizeGoal);

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
				System.out.println(getGoalNames());
			}
		};
		registerGoal(showGoals);

		Goal<GSProject> all = new Goal<GSProject>(project, "genall", kMerFastqTrieFileGoal, bloomFilterFileGoal) {
			@Override
			public boolean isMade() {
				return false;
			}

			@Override
			public void makeThis() {
			}
		};
		registerDefaultGoal(all);

		File fastq = project.getFastqFile();
		if (fastq != null) {
			Goal<GSProject> filterGoal = new FilterGoal(project, "filter", fastq,
					project.getConfig().isWriteDumpedFastq(), bloomFilterFileGoal, projectSetupGoal);
			registerGoal(filterGoal);

			Goal<GSProject> classifyGoal = new ClassifyGoal(project, "match", fastq, trieGoal, projectSetupGoal);
			registerGoal(classifyGoal);

			Goal<GSProject> krakenResCountGoal = new KrakenResCountGoal(project, "krakenres", fastq, taxNodesGoal,
					projectSetupGoal);
			registerGoal(krakenResCountGoal);

			Goal<GSProject> krakenResCountAllGoal = new KrakenResCountGoal(project, "krakenresall", fastq, null,
					projectSetupGoal);
			registerGoal(krakenResCountAllGoal);

			registerGoal(krakenResCountGoal);
			KrakenOutGoal fastqKrakenOutGoal = new KrakenOutGoal(project, "fastqkrakenout", fastq, projectSetupGoal);
			registerGoal(fastqKrakenOutGoal);

			ObjectGoal<KMerTrie<TaxidWithCount>, GSProject> trieFromKrakenResGoal = new TrieFromKrakenResGoal(project,
					"triefromkrakenres", fastq, taxNodesGoal, fastaFilesGoal, fastqKrakenOutGoal, fastaDownloadGoal,
					projectSetupGoal);
			registerGoal(trieFromKrakenResGoal);

			Goal<GSProject> krakenResErrorGoal = new KrakenResErrorGoal(project, "krakenerr", fastq,
					trieFromKrakenResGoal, projectSetupGoal);
			registerGoal(krakenResErrorGoal);
		}
	}
}