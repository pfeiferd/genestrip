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
package org.metagene.genestrip.goals.refseq;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.metagene.genestrip.ExecutionContext;
import org.metagene.genestrip.GSConfigKey;
import org.metagene.genestrip.GSGoalKey;
import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.refseq.AbstractRefSeqFastaReader;
import org.metagene.genestrip.refseq.AbstractStoreFastaReader;
import org.metagene.genestrip.refseq.AccessionMap;
import org.metagene.genestrip.refseq.RefSeqCategory;
import org.metagene.genestrip.tax.Rank;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;

/**
 * Goal that counts the total number of included k-mers (duplicates included) across the selected
 * RefSeq FASTA files, giving a rough estimate of the raw database size.
 *
 * @param <P> the project type
 */
public class FillSizeGoal<P extends GSProject> extends FastaReaderGoal<Long, P> {
	private final ObjectGoal<AccessionMap, P> accessionMapGoal;
	private final List<MyFastaReader> readers;

	/**
	 * Creates the goal, wiring the accession-map goal alongside the FASTA inputs it counts.
	 *
	 * @param project the project
	 * @param bundle the execution context
	 * @param categoriesGoal goal providing the selected RefSeq categories
	 * @param taxNodesGoal goal providing the tax id nodes to include
	 * @param fnaFilesGoal goal providing the downloaded RefSeq FASTA files
	 * @param additionalGoal goal providing additional FASTA files mapped to tax id nodes
	 * @param accessionMapGoal goal providing the accession-to-tax-id map
	 * @param deps additional goal dependencies
	 */
	@SafeVarargs
	public FillSizeGoal(P project, ExecutionContext bundle, ObjectGoal<Set<RefSeqCategory>, P> categoriesGoal,
						ObjectGoal<Set<TaxIdNode>, P> taxNodesGoal, RefSeqFnaFilesDownloadGoal fnaFilesGoal,
						ObjectGoal<Map<File, TaxIdNode>, P> additionalGoal,
						ObjectGoal<AccessionMap, P> accessionMapGoal, Goal<P>... deps) {
		super(project, GSGoalKey.FILLSIZE, bundle, categoriesGoal, taxNodesGoal, fnaFilesGoal, additionalGoal, Goal.append(deps, accessionMapGoal));
		this.accessionMapGoal = accessionMapGoal;
		readers = new ArrayList<>();
	}

	@Override
	protected void doMakeThis() {
		try {
			readFastas();
			long counter = 0;
			long dustSum = 0;
			long totalKmerSum = 0;

			for (MyFastaReader reader : readers) {
				counter += reader.getIncludedKmers();
				dustSum += reader.getDustCounter();
				totalKmerSum += reader.getTotalKmers();
			}
			set(counter);
			if (getLogger().isInfoEnabled()) {
				getLogger().info("All included kmers with duplicates: " + counter);
				getLogger().info("Estimated DB size in MB (without Bloom filter, with duplicates): " + (counter * 10) / (1024 * 1024) );
				if (intConfigValue(GSConfigKey.MAX_DUST) >= 0) {
					getLogger().info("Dust ratio: " + ((double) dustSum) / totalKmerSum);
				}
			}
		} catch (IOException e) {
			throw new RuntimeException(e);
		} finally {
			readers.clear();
			cleanUpThreads();
		}
	}

	@Override
	protected AbstractStoreFastaReader createFastaReader(AbstractRefSeqFastaReader.StringLong2DigitTrie regionsPerTaxid) {
		MyFastaReader fastaReader = new MyFastaReader(intConfigValue(GSConfigKey.FASTA_LINE_SIZE_BYTES),
				taxNodesGoal.get(), isIncludeRefSeqFna() ? accessionMapGoal.get() : null, intConfigValue(GSConfigKey.KMER_SIZE),
				intConfigValue(GSConfigKey.MAX_GENOMES_PER_TAXID),
				(Rank) configValue(GSConfigKey.MAX_GENOMES_PER_TAXID_RANK),
				longConfigValue(GSConfigKey.MAX_KMERS_PER_TAXID),
				intConfigValue(GSConfigKey.MAX_DUST),
				intConfigValue(GSConfigKey.STEP_SIZE),
				booleanConfigValue(GSConfigKey.COMPLETE_GENOMES_ONLY),
				regionsPerTaxid,
				booleanConfigValue(GSConfigKey.ENABLE_LOWERCASE_BASES));
		readers.add(fastaReader);
		return fastaReader;
	}

	/**
	 * FASTA reader that only counts k-mers (included, total and dust) without storing them.
	 */
	protected static class MyFastaReader extends AbstractStoreFastaReader {
		/**
		 * Creates a counting FASTA reader.
		 *
		 * @param bufferSize the read buffer size in bytes
		 * @param taxNodes the tax id nodes to include
		 * @param accessionMap the accession-to-tax-id map, or {@code null} if RefSeq FASTA is excluded
		 * @param k the k-mer length
		 * @param maxGenomesPerTaxId the maximum number of genomes per tax id
		 * @param maxGenomesPerTaxIdRank the rank at which the genome limit applies
		 * @param maxKmersPerTaxId the maximum number of k-mers per tax id
		 * @param maxDust the maximum dust value, or a negative value to disable dust filtering
		 * @param stepSize the k-mer step size
		 * @param completeGenomesOnly whether only complete-genome accessions are considered
		 * @param regionsPerTaxid the per-tax-id genome region trie
		 * @param enableLowerCaseBases whether lowercase bases are accepted
		 */
		public MyFastaReader(int bufferSize, Set<TaxIdNode> taxNodes, AccessionMap accessionMap, int k,
				int maxGenomesPerTaxId, Rank maxGenomesPerTaxIdRank, long maxKmersPerTaxId, int maxDust, int stepSize, boolean completeGenomesOnly, StringLong2DigitTrie regionsPerTaxid, boolean enableLowerCaseBases) {
			super(bufferSize, taxNodes, accessionMap, k, maxGenomesPerTaxId, maxGenomesPerTaxIdRank, maxKmersPerTaxId,
					maxDust, stepSize, completeGenomesOnly, regionsPerTaxid, enableLowerCaseBases);
		}

		/**
		 * Returns the number of included k-mers counted so far.
		 *
		 * @return the number of included k-mers counted so far
		 */
		public long getIncludedKmers() {
			return includedKmers;
		}

		/**
		 * Returns the total number of k-mers counted so far.
		 *
		 * @return the total number of k-mers counted so far
		 */
		public long getTotalKmers() {
			return totalKmers;
		}

		/**
		 * Returns the accumulated dust count.
		 *
		 * @return the accumulated dust count
		 */
		public long getDustCounter() {
			return dustCounter;
		}

		@Override
		protected boolean handleStore() {
			return true;
		}
	}
}