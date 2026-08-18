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
import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import net.agkn.hll.HLL;

import org.metagene.genestrip.ExecutionContext;
import org.metagene.genestrip.GSConfigKey;
import org.metagene.genestrip.GSGoalKey;
import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.refseq.AbstractRefSeqFastaReader;
import org.metagene.genestrip.util.MurmurHash3DropIn;
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
public class FillSizeGoal<P extends GSProject> extends FastaReaderGoal<FillSizeGoal.KMerCounts, P> {
	/**
	 * How many k-mers a fill would see, counted both ways: once as they come, duplicates and all, and
	 * once as distinct values.
	 * <p>
	 * The two differ by a lot - a k-mer shared by many genomes of the same species is counted once
	 * for every one of them in the first number and once altogether in the second - and which of them
	 * is wanted depends on the question. Sizing a structure that stores every k-mer separately calls
	 * for the first; sizing one that deduplicates, such as the temporary Bloom filter, calls for the
	 * second, and using the first there merely wastes memory.
	 */
	public static class KMerCounts implements Serializable {
		private static final long serialVersionUID = 1L;

		private final long withDuplicates;
		private final long distinct;

		/**
		 * Creates the pair of counts.
		 *
		 * @param withDuplicates the number of k-mers counted as they come, duplicates included
		 * @param distinct the estimated number of distinct k-mers among them
		 */
		public KMerCounts(long withDuplicates, long distinct) {
			this.withDuplicates = withDuplicates;
			this.distinct = distinct;
		}

		/**
		 * Returns the number of k-mers a fill would see, duplicates included. This one is exact.
		 *
		 * @return the number of k-mers including duplicates
		 */
		public long getWithDuplicates() {
			return withDuplicates;
		}

		/**
		 * Returns the estimated number of distinct k-mers among them.
		 * <p>
		 * Estimated, not counted: holding every k-mer seen would need the very memory this figure is
		 * meant to size. It comes from a HyperLogLog sketch of some 24 KB per reader thread, whose
		 * relative error is around half a percent and which may fall on either side of the truth.
		 *
		 * @return the estimated number of distinct k-mers
		 */
		public long getDistinct() {
			return distinct;
		}

		@Override
		public String toString() {
			return withDuplicates + " k-mers with duplicates, " + distinct + " distinct (estimated)";
		}
	}

	// The sketches of the reader threads are merged at the end, so every reader needs to hash a k-mer
	// to the same value: one base for all of them, fixed rather than random so that two runs over the
	// same data give the same size.
	private static final long HASH_BASE = 0x2545F4914F6CDD1DL;
	// Sizing of the HyperLogLog sketch: 2^15 registers of 6 bits, i.e. some 24 KB once it is fully
	// materialised, for a relative error of about half a percent. Six bits rather than the customary
	// five because a register has to count the leading zeroes of the largest run seen, and these
	// sketches take billions of k-mers.
	private static final int HLL_LOG2M = 15;
	private static final int HLL_REGISTER_WIDTH = 6;

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

			// Each reader sketched its own k-mers without any locking; merging the sketches afterwards
			// yields the same estimate as one sketch fed by all of them would have, because they hash
			// a k-mer alike and a register keeps a maximum, which does not care in what order or by
			// whom it was raised.
			HLL distinctSketch = newSketch();
			for (MyFastaReader reader : readers) {
				counter += reader.getIncludedKmers();
				dustSum += reader.getDustCounter();
				totalKmerSum += reader.getTotalKmers();
				distinctSketch.union(reader.getDistinctSketch());
			}
			long distinct = distinctSketch.cardinality();
			set(new KMerCounts(counter, distinct));
			if (getLogger().isInfoEnabled()) {
				getLogger().info("All included kmers with duplicates: " + counter);
				getLogger().info("Estimated distinct kmers: " + distinct);
				if (counter > 0) {
					getLogger().info("Duplication factor: " + ((double) counter) / distinct);
				}
				getLogger().info("Estimated DB size in MB (without Bloom filter, with duplicates): " + (counter * 10) / (1024 * 1024) );
				getLogger().info("Estimated DB size in MB (without Bloom filter, distinct only): " + (distinct * 10) / (1024 * 1024) );
				if (intConfigValue(GSConfigKey.MAX_DUST) >= 0) {
					getLogger().info("Dust ratio: " + ((double) dustSum) / totalKmerSum);
				}
			}
		} catch (IOException e) {
			throw new RuntimeException(e);
		} finally {
			readers.clear();
		}
	}

	/**
	 * Creates a HyperLogLog sketch of the sizing all readers share, so that theirs can be merged.
	 * <p>
	 * The sketch is left to promote itself through the library's representations - an exact list of
	 * values while there are few, a sparse map of registers next, the fully materialised registers in
	 * the end - rather than being forced to the last of them. That costs nothing here, since a sketch
	 * fed a database's worth of k-mers arrives at the full representation within its first moments,
	 * and it buys an <em>exact</em> count for a project small enough never to get there, where the
	 * full representation would have answered 101 to a hundred distinct k-mers.
	 *
	 * @return a new, empty sketch
	 */
	protected static HLL newSketch() {
		return new HLL(HLL_LOG2M, HLL_REGISTER_WIDTH);
	}

	@Override
	protected AbstractStoreFastaReader createFastaReader(AbstractRefSeqFastaReader.StringLong2DigitTrie regionsPerTaxid) {
		MyFastaReader fastaReader = new MyFastaReader(intConfigValue(GSConfigKey.FASTA_LINE_SIZE_BYTES),
				taxNodesGoal.get(), isIncludeRefSeqFna() ? accessionMapGoal.get() : null, intConfigValue(GSConfigKey.KMER_SIZE),
				intConfigValue(GSConfigKey.MAX_GENOMES_PER_TAXID),
				(Rank) configValue(GSConfigKey.MAX_GENOMES_PER_TAXID_RANK),
				longConfigValue(GSConfigKey.MAX_KMERS_PER_TAXID),
				intConfigValue(GSConfigKey.MAX_DUST),
				intConfigValue(GSConfigKey.KMER_SAMPLING),
				booleanConfigValue(GSConfigKey.ASSEMBLY_ACCESSIONS_ONLY),
				regionsPerTaxid,
				booleanConfigValue(GSConfigKey.ENABLE_LOWERCASE_BASES));
		readers.add(fastaReader);
		return fastaReader;
	}

	/**
	 * FASTA reader that only counts k-mers (included, total and dust) without storing them.
	 */
	protected static class MyFastaReader extends AbstractStoreFastaReader {
		// One sketch per reader, merged by the goal once the readers are done. A single shared sketch
		// would have to be locked for every k-mer, which is the hottest path there is here.
		private final HLL distinctSketch = newSketch();

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
		 * @param kMerSampling the k-mer step size
		 * @param assemblyAccessionsOnly whether only complete-genome accessions are considered
		 * @param regionsPerTaxid the per-tax-id genome region trie
		 * @param enableLowerCaseBases whether lowercase bases are accepted
		 */
		public MyFastaReader(int bufferSize, Set<TaxIdNode> taxNodes, AccessionMap accessionMap, int k,
				int maxGenomesPerTaxId, Rank maxGenomesPerTaxIdRank, long maxKmersPerTaxId, int maxDust, int kMerSampling, boolean assemblyAccessionsOnly, StringLong2DigitTrie regionsPerTaxid, boolean enableLowerCaseBases) {
			super(bufferSize, taxNodes, accessionMap, k, maxGenomesPerTaxId, maxGenomesPerTaxIdRank, maxKmersPerTaxId,
					maxDust, kMerSampling, assemblyAccessionsOnly, regionsPerTaxid, enableLowerCaseBases);
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

		/**
		 * Returns this reader's sketch of the distinct k-mers it saw, to be merged with the others.
		 *
		 * @return the sketch of this reader's k-mers
		 */
		public HLL getDistinctSketch() {
			return distinctSketch;
		}

		@Override
		protected boolean handleStore(long kmer) {
			// Hashing first matters: HyperLogLog reads its register index and its leading zeroes out of
			// the bits it is given, and a k-mer encodes two bits per base, so neighbouring k-mers would
			// otherwise land in neighbouring registers and the estimate would come out badly.
			distinctSketch.addRaw(MurmurHash3DropIn.hash64(kmer, HASH_BASE));
			return true;
		}
	}
}