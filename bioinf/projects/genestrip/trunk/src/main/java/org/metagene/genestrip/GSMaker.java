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
import java.io.PrintStream;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

import org.metagene.genestrip.bloom.FastqBloomFilter;
import org.metagene.genestrip.bloom.KMerBloomIndex;
import org.metagene.genestrip.goals.AssemblyFileDownloadGoal;
import org.metagene.genestrip.goals.BloomFilterFileGoal;
import org.metagene.genestrip.goals.BloomFilterSizeGoal;
import org.metagene.genestrip.goals.FastaFileDownloadGoal;
import org.metagene.genestrip.goals.FastasSizeGoal;
import org.metagene.genestrip.goals.KMerFastqGoal;
import org.metagene.genestrip.goals.KMerTrieFileGoal;
import org.metagene.genestrip.goals.KrakenFastqFileGoal;
import org.metagene.genestrip.goals.KrakenOutGoal;
import org.metagene.genestrip.goals.KrakenResCountGoal;
import org.metagene.genestrip.goals.TaxIdFileDownloadGoal;
import org.metagene.genestrip.goals.TrieFromKrakenResGoal;
import org.metagene.genestrip.goals.TrieFromKrakenResGoal.TaxidWithCount;
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
import org.metagene.genestrip.trie.KMerTrie.KMerTrieVisitor;
import org.metagene.genestrip.util.ByteArrayToString;
import org.metagene.genestrip.util.CountingDigitTrie;
import org.metagene.genestrip.util.StreamProvider;

public class GSMaker extends Maker<GSProject> {
	public GSMaker(GSProject project) {
		super(project);
	}

	protected void createGoals(GSProject project) {
		List<File> projectDirs = Arrays.asList(project.getFastasDir(), project.getFastqsDir(), project.getFiltersDir(),
				project.getKrakenOutDir(), project.getResultsDir());

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

		FileGoal<GSProject> assemblyGoal = new AssemblyFileDownloadGoal(project, "assemblydownload");
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
							getProject().getConfig().getCommonBaseDir(), getProject().getConfig().isUseGenBank(),
							taxTreeGoal.get());
					Map<TaxIdNode, List<FTPEntryWithQuality>> entries = assemblySummaryReader
							.getRelevantEntries(taxNodesGoal.get());
					set(entries);
				} catch (IOException e) {
					throw new RuntimeException(e);
				}
			}
		};

		FastaFileDownloadGoal fastaDownloadGoal = new FastaFileDownloadGoal(project, "fastasdownload", fastaFilesGoal,
				projectSetupGoal);
		registerGoal(fastaDownloadGoal);

		FastasSizeGoal fastasSizeGoal = new FastasSizeGoal(project, "fastassize", fastaDownloadGoal);

		KMerFastqGoal kmerFastqGoal = new KMerFastqGoal(project, "kmerfastqgen", fastasSizeGoal, fastaDownloadGoal,
				projectSetupGoal, fastaDownloadGoal);
		registerGoal(kmerFastqGoal);

		Goal<GSProject> krakenOutGoal = new KrakenOutGoal(project, "kmerkrakenout", project.getKmerFastqFile(),
				project.getKrakenOutFile(), kmerFastqGoal);
		registerGoal(krakenOutGoal);

		Goal<GSProject> trieGoal = new KMerTrieFileGoal(project, "triegen", taxNodesGoal, krakenOutGoal, kmerFastqGoal,
				projectSetupGoal);
		registerGoal(trieGoal);

		Goal<GSProject> krakenFastqGoal = new KrakenFastqFileGoal(project, "krakenfastq", taxNodesGoal,
				projectSetupGoal, taxNodesGoal, krakenOutGoal, kmerFastqGoal);
		registerGoal(krakenFastqGoal);

		BloomFilterSizeGoal bloomFilterSizeGoal = new BloomFilterSizeGoal(project, "bloomsize", taxNodesGoal,
				krakenOutGoal);
		registerGoal(bloomFilterSizeGoal);

		Goal<GSProject> bloomFilterFileGoal = new BloomFilterFileGoal(project, "bloomgen", bloomFilterSizeGoal,
				taxNodesGoal, krakenOutGoal, kmerFastqGoal, taxNodesGoal, bloomFilterSizeGoal);
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

			Goal<GSProject> classifyGoal = new FileListGoal<GSProject>(project, "classify",
					project.getTaxCountsFile("classify"), trieGoal) {
				@Override
				protected void makeFile(File file) {
					try {
						@SuppressWarnings("unchecked")
						KMerTrie<String> trie = (KMerTrie<String>) KMerTrie.load(getProject().getTrieFile());
						Map<String, Long> res = new FastqTrieClassifier(trie,
								getProject().getConfig().getMaxReadSizeBytes())
										.runClassifier(getProject().getFastqFile());
						PrintStream out = new PrintStream(StreamProvider.getOutputStreamForFile(file));
						CountingDigitTrie.print(res, out);
						out.close();
					} catch (IOException | ClassNotFoundException e) {
						throw new RuntimeException(e);
					}
				}
			};
			registerGoal(classifyGoal);

			Goal<GSProject> krakenResCountGoal = project.isUseKrakenOutFilter()
					? new KrakenResCountGoal(project, "krakenrescount", taxNodesGoal, taxNodesGoal, projectSetupGoal)
					: new KrakenResCountGoal(project, "krakenrescount", null, projectSetupGoal);
			registerGoal(krakenResCountGoal);

			Goal<GSProject> fastqKrakenOutGoal = new KrakenOutGoal(project, "fastqkrakenout", project.getFastqFile(),
					project.getFastqKrakenOutFile());
			registerGoal(fastqKrakenOutGoal);

			ObjectGoal<KMerTrie<TaxidWithCount>, GSProject> trieFromKrakenResGoal = new TrieFromKrakenResGoal(project, taxNodesGoal, fastaFilesGoal,
					taxNodesGoal, fastaFilesGoal, fastqKrakenOutGoal);
			registerGoal(trieFromKrakenResGoal);
			
			Goal<GSProject> krakenResErrorGoal = new FileListGoal<GSProject>(project, "krakenreserr", project.getKrakenErrFile(), trieFromKrakenResGoal) {
				@Override
				protected void makeFile(File file) throws IOException {
					final PrintStream ps = new PrintStream(StreamProvider.getOutputStreamForFile(file));
					
					ps.println("taxid;count;kmer");
					KMerTrie<TaxidWithCount> trie = trieFromKrakenResGoal.get();
					Map<String, Integer> errPerTaxid = new HashMap<String, Integer>();
					
					trie.visit(new KMerTrieVisitor<TrieFromKrakenResGoal.TaxidWithCount>() {						
						@Override
						public void nextValue(KMerTrie<TaxidWithCount> trie, byte[] kmer, TaxidWithCount value) {
							if (value != null) {
								ps.print(value.getTaxid());
								ps.print(';');
								ps.print(value.getCount());
								ps.print(';');
								ByteArrayToString.print(kmer, ps);
								ps.println(';');
								Integer res = errPerTaxid.get(value.getTaxid());
								res = (res == null) ? value.getCount() : res + value.getCount();
								errPerTaxid.put(value.getTaxid(), res);
							}
						}
					}, false);
					ps.println("taxid;aggregated count");
					for (Entry<String, Integer> entry : errPerTaxid.entrySet()) {
						ps.print(entry.getKey());
						ps.print(';');
						ps.println(entry.getValue());
						ps.println(';');
					}
					
					ps.close();
				}
			};
			registerGoal(krakenResErrorGoal);
		}
	}
}