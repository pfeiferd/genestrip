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

public class FillSizeGoal extends FastaReaderGoal<Long> {
	private final ObjectGoal<AccessionMap, GSProject> accessionMapGoal;

	private final List<MyFastaReader> readers;

	@SafeVarargs
	public FillSizeGoal(GSProject project, ExecutionContext bundle, ObjectGoal<Set<RefSeqCategory>, GSProject> categoriesGoal,
						ObjectGoal<Set<TaxIdNode>, GSProject> taxNodesGoal, RefSeqFnaFilesDownloadGoal fnaFilesGoal,
						ObjectGoal<Map<File, TaxIdNode>, GSProject> additionalGoal,
						ObjectGoal<AccessionMap, GSProject> accessionMapGoal, Goal<GSProject>... deps) {
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
				regionsPerTaxid);
		readers.add(fastaReader);
		return fastaReader;
	}

	protected static class MyFastaReader extends AbstractStoreFastaReader {
		public MyFastaReader(int bufferSize, Set<TaxIdNode> taxNodes, AccessionMap accessionMap, int k,
				int maxGenomesPerTaxId, Rank maxGenomesPerTaxIdRank, long maxKmersPerTaxId, int maxDust, int stepSize, boolean completeGenomesOnly, StringLong2DigitTrie regionsPerTaxid) {
			super(bufferSize, taxNodes, accessionMap, k, maxGenomesPerTaxId, maxGenomesPerTaxIdRank, maxKmersPerTaxId,
					maxDust, stepSize, completeGenomesOnly, regionsPerTaxid);
		}

		public long getIncludedKmers() {
			return includedKmers;
		}

		public long getTotalKmers() {
			return totalKmers;
		}

		public long getDustCounter() {
			return dustCounter;
		}

		@Override
		protected boolean handleStore() {
			return true;
		}

	/* This would be faster but the estimate in kmers is a little less accurate:
	@Override
	protected void dataLine() {
		if (includeRegion) {
			if (isAllowMoreKmers()) {
				bpsInRegion += size - 1;
				if (bpsInRegion >= k) {
					kmersInRegion = ((bpsInRegion - k + 1) / stepSize); // TODO not yet correct
				}
			}
		}
	}
	 */
	}
}