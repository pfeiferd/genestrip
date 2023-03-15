package org.metagene.genestrip.gen;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.Set;

import org.metagene.genestrip.bloom.FastqBloomFilter;
import org.metagene.genestrip.bloom.KMerBloomIndex;
import org.metagene.genestrip.fastqgen.KrakenExecutor;
import org.metagene.genestrip.gen.goals.AssemblyFileDownloadGoal;
import org.metagene.genestrip.gen.goals.BloomFilterFileGoal;
import org.metagene.genestrip.gen.goals.FastaFileDownloadGoal;
import org.metagene.genestrip.gen.goals.FastasSizeGoal;
import org.metagene.genestrip.gen.goals.KMerFastqGoal;
import org.metagene.genestrip.gen.goals.KMerTrieFileGoal;
import org.metagene.genestrip.gen.goals.KrakenFastqFileGoal;
import org.metagene.genestrip.gen.goals.TaxIdFileDownloadGoal;
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

public class Generator extends Maker<Project> {
	public Generator(Project project) {
		super(project);
	}

	protected void createGoals(Project project) {
		Config config = project.getConfig();

		List<File> projectDirs = Arrays.asList(project.getFastasDir(), project.getFastqsDir(), project.getFiltersDir(),
				project.getResultsDir());

		Goal<Project> projectSetupGoal = new FileListGoal<Project>(project, "setup", projectDirs) {
			@Override
			protected void makeFile(File file) throws IOException {
				file.mkdir();
			}
		};
		registerGoal(projectSetupGoal);

		Goal<Project> taxDBGoal = new TaxIdFileDownloadGoal(project, "taxdownload");
		registerGoal(taxDBGoal);

		ObjectGoal<TaxTree, Project> taxTreeGoal = new ObjectGoal<TaxTree, Project>(project, "taxtree", taxDBGoal) {
			@Override
			public void makeThis() {
				set(new TaxTree(config.getCommonBaseDir()));
			}
		};
		registerGoal(taxTreeGoal);

		ObjectGoal<Set<TaxIdNode>, Project> taxNodesGoal = new ObjectGoal<Set<TaxIdNode>, Project>(project, "taxids",
				taxTreeGoal) {
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
					}
					set(taxIdNodes);
				} catch (IOException e) {
					throw new RuntimeException(e);
				}
			}
		};
		registerGoal(taxNodesGoal);

		FileGoal<Project> assemblyGoal = new AssemblyFileDownloadGoal(project);
		registerGoal(assemblyGoal);

		Goal<Project> commonGoal = new Goal<Project>(project, "common", assemblyGoal, taxDBGoal) {
			@Override
			public boolean isMade() {
				return false;
			}
		};
		registerGoal(commonGoal);

		ObjectGoal<List<FTPEntryWithQuality>, Project> fastaFilesGoal = new ObjectGoal<List<FTPEntryWithQuality>, Project>(
				project, "fastafiles", assemblyGoal, taxNodesGoal) {
			@Override
			public void makeThis() {
				try {
					AssemblySummaryReader assemblySummaryReader = new AssemblySummaryReader(config.getCommonBaseDir(),
							taxTreeGoal.get());
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

		Goal<Project> krakenOutGoal = new FileListGoal<Project>(project, "krakenout", kmerFastqGoal, projectSetupGoal) {
			@Override
			protected void makeFile(File krakenOut) {
				File fastq = getProject().getKmerFastqFile();
				KrakenExecutor krakenExecutor = new KrakenExecutor(config.getKrakenBinFolder(),
						config.getKrakenExecExpr());
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

		Goal<Project> trieGoal = new KMerTrieFileGoal(project, taxNodesGoal, krakenOutGoal, kmerFastqGoal,
				projectSetupGoal);
		registerGoal(trieGoal);

		Goal<Project> krakenFastqGoal = new KrakenFastqFileGoal(project, taxNodesGoal, projectSetupGoal);
		registerGoal(krakenFastqGoal);

		Goal<Project> bloomFilterFileGoal = new BloomFilterFileGoal(project, kmerFastqGoal, taxNodesGoal,
				projectSetupGoal, kmerFastqGoal);
		registerGoal(bloomFilterFileGoal);

		Goal<Project> showGoals = new Goal<Project>(project, "show") {
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

		Goal<Project> all = new Goal<Project>(project, "genall", trieGoal, krakenFastqGoal, bloomFilterFileGoal) {
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
			Goal<Project> filterGoal = new FileListGoal<Project>(project, "filter", project.getFilteredFastqFile(),
					bloomFilterFileGoal) {
				@Override
				protected void makeFile(File file) {
					try {
						new FastqBloomFilter(KMerBloomIndex.load(getProject().getBloomFilterFile()),
								config.getPosRatioFilter(), config.getMinPosCountFilter(), config.getMaxReadSizeBytes())
										.runFilter(getProject().getFastqFile(), file, null);
					} catch (ClassNotFoundException e) {
						throw new RuntimeException(e);
					} catch (IOException e) {
						throw new RuntimeException(e);
					}
				}
			};
			registerGoal(filterGoal);
		}
	}
}