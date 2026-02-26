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
import java.util.*;

import net.agkn.hll.HLL;
import net.agkn.hll.HLLType;
import org.apache.commons.codec.language.MatchRatingApproachEncoder;
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
import org.metagene.genestrip.util.MurmurHash3DropIn;

import static org.apache.commons.math3.util.MathUtils.hash;

public class FillSizeGoal<P extends GSProject> extends FastaReaderGoal<Long, P> {
	private final ObjectGoal<AccessionMap, P> accessionMapGoal;
	private final List<MyFastaReader> readers;

	private Random random;
	private HLL hll;
	private final boolean multiThreading;

	@SafeVarargs
	public FillSizeGoal(P project, ExecutionContext bundle, ObjectGoal<Set<RefSeqCategory>, P> categoriesGoal,
						ObjectGoal<Set<TaxIdNode>, P> taxNodesGoal, RefSeqFnaFilesDownloadGoal fnaFilesGoal,
						ObjectGoal<Map<File, TaxIdNode>, P> additionalGoal,
						ObjectGoal<AccessionMap, P> accessionMapGoal, Goal<P>... deps) {
		super(project, GSGoalKey.FILLSIZE, bundle, categoriesGoal, taxNodesGoal, fnaFilesGoal, additionalGoal, Goal.append(deps, accessionMapGoal));
		this.accessionMapGoal = accessionMapGoal;
		readers = new ArrayList<>();
		multiThreading = bundle.getThreads() > 0;
		random = new Random(42);
	}

	@Override
	protected void doMakeThis() {
		try {
			hll = new HLL(15, 6, -1, false, HLLType.FULL);
			readFastas();
			long counter = 0;
			long dustSum = 0;
			long totalKmerSum = 0;

			set(hll.cardinality());
			if (getLogger().isInfoEnabled()) {
				for (MyFastaReader reader : readers) {
					counter += reader.getIncludedKmers();
					dustSum += reader.getDustCounter();
					totalKmerSum += reader.getTotalKmers();
				}
				long dedup = get();
				getLogger().info("All included kmers with duplicates: " + counter);
				getLogger().info("Deduplicated kmers: " + dedup);
				getLogger().info("Duplication factor: " + ((double) counter) / dedup);
				getLogger().info("Estimated DB size in MB (without Bloom filter): " + (dedup * 10) / (1024 * 1024) );
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
				regionsPerTaxid,
				booleanConfigValue(GSConfigKey.ENABLE_LOWERCASE_BASES));
		readers.add(fastaReader);
		return fastaReader;
	}

	protected class MyFastaReader extends AbstractStoreFastaReader {
		private final long hashBase;

		public MyFastaReader(int bufferSize, Set<TaxIdNode> taxNodes, AccessionMap accessionMap, int k,
				int maxGenomesPerTaxId, Rank maxGenomesPerTaxIdRank, long maxKmersPerTaxId, int maxDust, int stepSize, boolean completeGenomesOnly, StringLong2DigitTrie regionsPerTaxid, boolean enableLowerCaseBases) {
			super(bufferSize, taxNodes, accessionMap, k, maxGenomesPerTaxId, maxGenomesPerTaxIdRank, maxKmersPerTaxId,
					maxDust, stepSize, completeGenomesOnly, regionsPerTaxid, enableLowerCaseBases);
			hashBase = random.nextLong();
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
			// Hashing proves to be extremely important for a good estimate.
			long hash = MurmurHash3DropIn.hash64(byteRingBuffer.getStandardKMer(), hashBase);
			if (multiThreading) {
				synchronized (hll) {
					hll.addRaw(hash);
				}
			}
			else {
				hll.addRaw(hash);
			}
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