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
import java.util.Map;
import java.util.Set;

import org.metagene.genestrip.ExecutionContext;
import org.metagene.genestrip.GSConfigKey;
import org.metagene.genestrip.GSGoalKey;
import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.bloom.MurmurCGATBloomFilter;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.refseq.AbstractRefSeqFastaReader;
import org.metagene.genestrip.refseq.AbstractStoreFastaReader;
import org.metagene.genestrip.refseq.AccessionMap;
import org.metagene.genestrip.refseq.RefSeqCategory;
import org.metagene.genestrip.tax.Rank;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;

public class FillBloomFilterGoal extends FastaReaderGoal<MurmurCGATBloomFilter> {
	private final ObjectGoal<AccessionMap, GSProject> accessionMapGoal;
	private final ObjectGoal<Long, GSProject> sizeGoal;

	private MurmurCGATBloomFilter filter;

	private final boolean multiThreading;

	@SafeVarargs
	public FillBloomFilterGoal(GSProject project, ExecutionContext bundle, ObjectGoal<Set<RefSeqCategory>, GSProject> categoriesGoal,
							   ObjectGoal<Set<TaxIdNode>, GSProject> taxNodesGoal, RefSeqFnaFilesDownloadGoal fnaFilesGoal,
							   ObjectGoal<Map<File, TaxIdNode>, GSProject> additionalGoal,
							   ObjectGoal<AccessionMap, GSProject> accessionMapGoal, FillSizeGoal sizeGoal, Goal<GSProject>... deps) {
		super(project, GSGoalKey.TEMPINDEX, bundle, categoriesGoal, taxNodesGoal, fnaFilesGoal, additionalGoal, Goal.append(deps, sizeGoal));
		this.accessionMapGoal = accessionMapGoal;
		this.sizeGoal = sizeGoal;
		multiThreading = bundle.getThreads() > 0;
	}
	
	@Override
	protected void allDependentsMade() {
		// To save memory...
		doCleanThis();
	}

	@Override
	protected void doMakeThis() {
		try {
			filter = new MurmurCGATBloomFilter(intConfigValue(GSConfigKey.KMER_SIZE),
					doubleConfigValue(GSConfigKey.TEMP_BLOOM_FILTER_FPP));
			filter.ensureExpectedSize(sizeGoal.get(), false);
			readFastas();
			set(filter);
			if (getLogger().isWarnEnabled()) {
				getLogger().warn("Final Bloom filter and store size in kmers: " + filter.getEntries());
				getLogger().warn("Approx. DB Size in MB (without Bloom filter): " + (filter.getEntries() * 10) / (1024 * 1024));
			}
		} catch (IOException e) {
			throw new RuntimeException(e);
		} finally {
			filter = null;
			cleanUpThreads();
		}
	}

	@Override
	protected AbstractRefSeqFastaReader createFastaReader() {
		return new MyFastaReader(intConfigValue(GSConfigKey.FASTA_LINE_SIZE_BYTES),
				taxNodesGoal.get(),
				booleanConfigValue(GSConfigKey.REF_SEQ_DB) ? accessionMapGoal.get() : null, filter,
				intConfigValue(GSConfigKey.MAX_GENOMES_PER_TAXID),
				(Rank) configValue(GSConfigKey.MAX_GENOMES_PER_TAXID_RANK),
				longConfigValue(GSConfigKey.MAX_KMERS_PER_TAXID),
				intConfigValue(GSConfigKey.MAX_DUST),
				intConfigValue(GSConfigKey.STEP_SIZE),
				booleanConfigValue(GSConfigKey.COMPLETE_GENOMES_ONLY));
	}

	protected class MyFastaReader extends AbstractStoreFastaReader {
		private final MurmurCGATBloomFilter filter;

		public MyFastaReader(int bufferSize, Set<TaxIdNode> taxNodes, AccessionMap accessionMap,
							 MurmurCGATBloomFilter filter, int maxGenomesPerTaxId, Rank maxGenomesPerTaxIdRank, long maxKmersPerTaxId, int maxDust, int stepSize, boolean completeGenomesOnly) {
			super(bufferSize, taxNodes, accessionMap, filter.getK(), maxGenomesPerTaxId, maxGenomesPerTaxIdRank, maxKmersPerTaxId, maxDust, stepSize, completeGenomesOnly);
			this.filter = filter;
		}

		@Override
		protected boolean handleStore() {
			long kmer = byteRingBuffer.getStandardKMer();
			if (!filter.containsLong(kmer)) {
				if (multiThreading) {
					synchronized (filter) {
						// This is a trick to enable more parallelism -
						// check again after synchronized to avoid synchronized further outside...
						if (!filter.containsLong(kmer)) {
							filter.putLong(kmer);
							return true;
						}
					}
				}
				else {
					filter.putLong(kmer);
				}
			}
			return false;
		}

		@Override
		protected void done() {
			super.done();
			if (getLogger().isInfoEnabled()) {
				long entries = filter.getEntries();
				getLogger().info("Number of included regions: " + includedCounter);
				getLogger().info("Total Bloom filter entries: " + entries);
				getLogger().info("Resulting approx. DB Size in MB (without Bloom filter): " + (entries * 10) / (1024 * 1024) );
			}
		}
	}
}