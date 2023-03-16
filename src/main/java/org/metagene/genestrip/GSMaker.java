package org.metagene.genestrip;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.metagene.genestrip.bloom.FastqBloomFilter;
import org.metagene.genestrip.bloom.KMerBloomIndex;
import org.metagene.genestrip.fastqgen.KrakenExecutor;
import org.metagene.genestrip.goals.AssemblyFileDownloadGoal;
import org.metagene.genestrip.goals.BloomFilterFileGoal;
import org.metagene.genestrip.goals.FastaFileDownloadGoal;
import org.metagene.genestrip.goals.FastasSizeGoal;
import org.metagene.genestrip.goals.KMerFastqGoal;
import org.metagene.genestrip.goals.KMerTrieFileGoal;
import org.metagene.genestrip.goals.KrakenFastqFileGoal;
import org.metagene.genestrip.goals.TaxIdFileDownloadGoal;
import org.metagene.genestrip.make.FileGoal;
import org.metagene.genestrip.make.FileListGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.Maker;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.tax.AssemblySummaryReader;
import org.metagene.genestrip.tax.AssemblySummaryReader.FTPEntryWithQuality;
import org.metagene.genestrip.tax.TaxIdCollector;
import org.metagene.genestrip.tax.TaxTree;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;
import org.metagene.genestrip.trie.FastqTrieClassifier;
import org.metagene.genestrip.trie.KMerTrie;
import org.metagene.genestrip.util.CountingDigitTrie;

public class GSMaker extends Maker<GSProject> {
	public GSMaker(GSProject project) {
		super(project);
	}

	protected void createGoals(GSProject project) {
		List<File> projectDirs = Arrays.asList(project.getFastasDir(), project.getFastqsDir(), project.getFiltersDir(),
				project.getResultsDir());

		Goal<GSProject> projectSetupGoal = new FileListGoal<GSProject>(project, "setup", projectDirs) {
			@Override
			protected void makeFile(File file) throws IOException {
				file.mkdir();
			}
		};
		registerGoal(projectSetupGoal);

		Goal<GSProject> taxDBGoal = new TaxIdFileDownloadGoal(project, "taxdownload");
		registerGoal(taxDBGoal);

		ObjectGoal<TaxTree, GSProject> taxTreeGoal = new ObjectGoal<TaxTree, GSProject>(project, "taxtree", taxDBGoal) {
			@Override
			public void makeThis() {
				set(new TaxTree(getProject().getConfig().getCommonBaseDir()));
			}
		};
		registerGoal(taxTreeGoal);

		ObjectGoal<Set<TaxIdNode>, GSProject> taxNodesGoal = new ObjectGoal<Set<TaxIdNode>, GSProject>(project,
				"taxids", taxTreeGoal) {
			@Override
			public void makeThis() {
				try {
					TaxIdCollector taxIdCollector = new TaxIdCollector(taxTreeGoal.get());
					Set<TaxIdNode> taxIdNodes = taxIdCollector.readFromFile(getProject().getTaxIdsFile());
					if (getLogger().isInfoEnabled()) {
						getLogger().info("Requested tax ids: " + taxIdNodes);
					}
					taxIdNodes = taxIdCollector.withDescendants(taxIdNodes);
					if (getLogger().isInfoEnabled()) {
						getLogger().info("Completed tax ids: " + taxIdNodes);
						getLogger().info("Number of completed tax ids: " + taxIdNodes.size());
					}
					set(taxIdNodes);
				} catch (IOException e) {
					throw new RuntimeException(e);
				}
			}
		};
		registerGoal(taxNodesGoal);

		FileGoal<GSProject> assemblyGoal = new AssemblyFileDownloadGoal(project);
		registerGoal(assemblyGoal);

		Goal<GSProject> commonGoal = new Goal<GSProject>(project, "common", assemblyGoal, taxDBGoal) {
			@Override
			public boolean isMade() {
				return false;
			}
		};
		registerGoal(commonGoal);

		ObjectGoal<List<FTPEntryWithQuality>, GSProject> fastaFilesGoal = new ObjectGoal<List<FTPEntryWithQuality>, GSProject>(
				project, "fastafiles", assemblyGoal, taxNodesGoal) {
			@Override
			public void makeThis() {
				try {
					AssemblySummaryReader assemblySummaryReader = new AssemblySummaryReader(
							getProject().getConfig().getCommonBaseDir(), taxTreeGoal.get());
					List<FTPEntryWithQuality> entries = assemblySummaryReader
							.getRelevantEntriesAsList(project.getFastaQuality(), taxNodesGoal.get());
					set(entries);
				} catch (IOException e) {
					throw new RuntimeException(e);
				}
			}
		};

		FastaFileDownloadGoal fastaDownloadGoal = new FastaFileDownloadGoal(project, fastaFilesGoal, projectSetupGoal);
		registerGoal(fastaDownloadGoal);

		FastasSizeGoal fastasSizeGoal = new FastasSizeGoal(project, fastaDownloadGoal);

		KMerFastqGoal kmerFastqGoal = new KMerFastqGoal(project, fastasSizeGoal, fastaDownloadGoal, projectSetupGoal,
				fastaDownloadGoal);
		registerGoal(kmerFastqGoal);

		Goal<GSProject> krakenOutGoal = new FileListGoal<GSProject>(project, "krakenout", project.getKrakenOutFile(),
				kmerFastqGoal, projectSetupGoal) {
			@Override
			protected void makeFile(File krakenOut) {
				File fastq = getProject().getKmerFastqFile();
				KrakenExecutor krakenExecutor = new KrakenExecutor(getProject().getConfig().getKrakenBinFolder(),
						getProject().getConfig().getKrakenExecExpr());
				if (getLogger().isInfoEnabled()) {
					String execLine = krakenExecutor.genExecLine(getProject().getKrakenDB(), fastq);
					getLogger().info("Run kraken with " + execLine);
				}
				try {
					krakenExecutor.execute(getProject().getKrakenDB(), fastq, krakenOut);
				} catch (InterruptedException | IOException e) {
					throw new RuntimeException(e);
				}
				if (getLogger().isInfoEnabled()) {
					getLogger().info("Finished kraken");
				}
			}
		};
		registerGoal(krakenOutGoal);

		Goal<GSProject> trieGoal = new KMerTrieFileGoal(project, taxNodesGoal, krakenOutGoal, kmerFastqGoal,
				projectSetupGoal);
		registerGoal(trieGoal);

		Goal<GSProject> krakenFastqGoal = new KrakenFastqFileGoal(project, taxNodesGoal, projectSetupGoal);
		registerGoal(krakenFastqGoal);

		Goal<GSProject> bloomFilterFileGoal = new BloomFilterFileGoal(project, kmerFastqGoal, taxNodesGoal,
				projectSetupGoal, kmerFastqGoal);
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

		Goal<GSProject> all = new Goal<GSProject>(project, "genall", trieGoal, krakenFastqGoal, bloomFilterFileGoal) {
			@Override
			public boolean isMade() {
				return false;
			}

			@Override
			public void makeThis() {
			}
		};
		registerDefaultGoal(all);

		if (project.getFastqFile() != null) {
			Goal<GSProject> filterGoal = new FileListGoal<GSProject>(project, "filter", project.getFilteredFastqFile(),
					bloomFilterFileGoal) {
				@Override
				protected void makeFile(File file) {
					try {
						new FastqBloomFilter(KMerBloomIndex.load(getProject().getBloomFilterFile()),
								getProject().getConfig().getPosRatioFilter(),
								getProject().getConfig().getMinPosCountFilter(),
								getProject().getConfig().getMaxReadSizeBytes()).runFilter(getProject().getFastqFile(),
										file, getProject().getDumpFastqFile());
					} catch (IOException | ClassNotFoundException e) {
						throw new RuntimeException(e);
					}
				}
			};
			registerGoal(filterGoal);

			Goal<GSProject> classifyGoal = new FileListGoal<GSProject>(project, "classify", project.getTaxCountsFile(),
					trieGoal) {
				@Override
				protected void makeFile(File file) {
					try {
						@SuppressWarnings("unchecked")
						KMerTrie<String> trie = (KMerTrie<String>) KMerTrie.load(getProject().getTrieFile());
						Map<String, Long> res = new FastqTrieClassifier(trie,
								getProject().getConfig().getMaxReadSizeBytes()).runClassifier(project.getFastqFile());
						FileOutputStream out = new FileOutputStream(project.getTaxCountsFile());
						CountingDigitTrie.print(res, out);
						out.close();
					} catch (IOException | ClassNotFoundException e) {
						throw new RuntimeException(e);
					}
				}
			};
			registerGoal(classifyGoal);
		}
	}
}