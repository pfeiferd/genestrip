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
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

import org.metagene.genestrip.ExecutionContext;
import org.metagene.genestrip.GSConfigKey;
import org.metagene.genestrip.GSGoalKey;
import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.GoalKey;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.refseq.AbstractRefSeqFastaReader;
import org.metagene.genestrip.refseq.AbstractStoreFastaReader;
import org.metagene.genestrip.refseq.AccessionMap;
import org.metagene.genestrip.refseq.RefSeqCategory;
import org.metagene.genestrip.store.Database;
import org.metagene.genestrip.store.KMerStore;
import org.metagene.genestrip.tax.Rank;
import org.metagene.genestrip.tax.SmallTaxTree;
import org.metagene.genestrip.tax.TaxTree;
import org.metagene.genestrip.tax.SmallTaxTree.SmallTaxIdNode;
import org.metagene.genestrip.tax.TaxTree;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;

/**
 * Goal that builds a per-species statistic of the taxonomic rank at which each of a species' k-mers
 * ended up in a Genestrip database. It re-reads the RefSeq (and additional) genomic FASTA files, and
 * for every k-mer of a species it looks that k-mer up in the loaded database's {@link KMerStore},
 * resolves the stored tax id to its {@link Rank} and tallies it.
 * <p>
 * The result is a {@code Map<String, long[]>} keyed by the species' tax id. Each value is an array
 * of counts indexed by {@link Rank#ordinal()}: entry {@code counts[rank.ordinal()]} holds the number
 * of the species' k-mers that landed on a node of that rank in the database. The array has
 * {@link Rank#values()}{@code .length} entries so that every rank has a fixed slot.
 * <p>
 * The species of a k-mer is determined by walking up from the region's tax node to its nearest
 * ancestor of rank {@link Rank#SPECIES}; if there is none, the region's own tax node is used. K-mers
 * not present in the database are ignored, as they did not land anywhere. Each qualifying k-mer
 * occurrence is counted (there is no de-duplication of k-mers that recur within a species' genomes).
 * <p>
 * A k-mer may land on a node whose rank has no well-defined level ({@link Rank#isIndeterminate()},
 * e.g. {@code clade} or {@code no rank}). Such a node cannot be placed on the rank axis by its rank
 * alone, so it is counted under the rank of its nearest ranked ancestor in the taxonomy tree (a
 * {@code no rank} node whose nearest ranked ancestor is a {@code genus}, for instance, is counted as
 * {@link Rank#GENUS}). This resolution is precomputed once from the database's taxonomy tree before
 * the multi-threaded run, so the count arrays only ever hold well-defined ranks and downstream binning
 * by rank works directly. A node with no ranked ancestor at all is not counted.
 *
 * @param <P> the project type
 */
public class KMerRankStatsGoal<P extends GSProject> extends FastaReaderGoal<Map<String, long[]>, P>
		implements Goal.LogHeapInfo {
	// Fixed length of every per-species count array: one slot per taxonomic rank, addressed by ordinal.
	private static final int RANK_COUNT = Rank.values().length;

	private final ObjectGoal<AccessionMap, P> accessionMapGoal;
	private final ObjectGoal<Database, P> dbGoal;

	// Populated for the duration of doMakeThis() and read by the (multi-threaded) fasta readers.
	private KMerStore<String> store;
	private SmallTaxTree taxTree;
	private Map<String, long[]> map;
	// For nodes whose own rank is indeterminate: the ordinal of the representative rank they are
	// counted under, resolved from their position in the tree. Precomputed before the parallel run so
	// the reader threads only read it. Determinate-rank nodes are absent (they use their own ordinal).
	private Map<String, Integer> effectiveOrdinalByTaxid;

	/**
	 * Creates the goal, wiring the accession map and the loaded database it looks k-mers up in, in
	 * addition to the standard FASTA-reader dependencies.
	 *
	 * @param project          the project this goal belongs to
	 * @param key              the key identifying this goal
	 * @param bundle           the execution context supplying the worker threads
	 * @param categoriesGoal   the goal supplying the RefSeq categories to read
	 * @param taxNodesGoal     the goal supplying the required taxonomy nodes
	 * @param fnaFilesGoal     the goal supplying the downloaded RefSeq {@code .fna} files
	 * @param additionalGoal   the goal supplying additional FASTA files mapped to their tax node
	 * @param accessionMapGoal the goal supplying the accession-to-taxid map
	 * @param dbGoal           the goal supplying the loaded database to look k-mers up in
	 * @param taxTreeGoal      the goal supplying the taxonomy tree; kept as a direct dependency so that
	 *                         the (aggressively cleaned) tree is not freed and rebuilt as a second
	 *                         instance while this goal runs - which would make the region nodes resolved
	 *                         via the accession map and the requested nodes from {@code taxNodesGoal}
	 *                         come from different tree instances, breaking their identity comparison
	 * @param deps             any further goals this goal depends on
	 */
	@SafeVarargs
	public KMerRankStatsGoal(P project, GoalKey key, ExecutionContext bundle,
							 ObjectGoal<Set<RefSeqCategory>, P> categoriesGoal,
							 ObjectGoal<Set<TaxIdNode>, P> taxNodesGoal, RefSeqFnaFilesDownloadGoal fnaFilesGoal,
							 ObjectGoal<Map<File, TaxIdNode>, P> additionalGoal,
							 ObjectGoal<AccessionMap, P> accessionMapGoal, ObjectGoal<Database, P> dbGoal,
							 ObjectGoal<TaxTree, P> taxTreeGoal, Goal<P>... deps) {
		super(project, key, bundle, categoriesGoal, taxNodesGoal, fnaFilesGoal, additionalGoal,
				Goal.append(deps, accessionMapGoal, dbGoal, taxTreeGoal));
		this.accessionMapGoal = accessionMapGoal;
		this.dbGoal = dbGoal;
	}

	@Override
	protected void doMakeThis() {
		try {
			map = new HashMap<>();
			effectiveOrdinalByTaxid = new HashMap<>();
			Database database = dbGoal.get();
			store = database.getKmerStore();
			taxTree = database.getTaxTree();
			// Two things are precomputed here, from the database's own taxonomy tree, before the
			// multi-threaded readFastas() run so that both maps are read-only (and thus safe to read
			// concurrently) while the reader threads work:
			// 1. A map entry for every species that can occur as a key, so the map is never structurally
			//    modified during the run - the readers only mutate an already-present long[] value, each
			//    guarded by its species node's own monitor (fine-grained per-species locking).
			// 2. For every node whose own rank is indeterminate (e.g. clade / no rank), the ordinal of
			//    the representative rank it is counted under, resolved from its position in the tree.
			for (SmallTaxIdNode node : taxTree) {
				map.computeIfAbsent(toSpeciesNode(node).getTaxId(), k -> new long[RANK_COUNT]);
				Rank rank = node.getRank();
				if (rank == null || rank.isIndeterminate()) {
					Rank effectiveRank = resolveEffectiveRank(node);
					if (effectiveRank != null) {
						effectiveOrdinalByTaxid.put(node.getTaxId(), effectiveRank.ordinal());
					}
				}
			}
			readFastas();
			set(map);
		} catch (IOException e) {
			throw new RuntimeException(e);
		} finally {
			map = null;
			store = null;
			taxTree = null;
			effectiveOrdinalByTaxid = null;
			cleanUpThreads();
		}
	}

	@Override
	protected AbstractStoreFastaReader createFastaReader(AbstractRefSeqFastaReader.StringLong2DigitTrie regionsPerTaxid) {
		return new MyFastaReader(intConfigValue(GSConfigKey.FASTA_LINE_SIZE_BYTES), taxNodesGoal.get(),
				isIncludeRefSeqFna() ? accessionMapGoal.get() : null, intConfigValue(GSConfigKey.KMER_SIZE),
				intConfigValue(GSConfigKey.MAX_GENOMES_PER_TAXID),
				(Rank) configValue(GSConfigKey.MAX_GENOMES_PER_TAXID_RANK),
				longConfigValue(GSConfigKey.MAX_KMERS_PER_TAXID), intConfigValue(GSConfigKey.MAX_DUST),
				intConfigValue(GSConfigKey.STEP_SIZE), booleanConfigValue(GSConfigKey.COMPLETE_GENOMES_ONLY),
				regionsPerTaxid, booleanConfigValue(GSConfigKey.ENABLE_LOWERCASE_BASES));
	}

	/**
	 * Walks up the taxonomy from the given node to the nearest ancestor of rank {@link Rank#SPECIES}.
	 *
	 * @param node the node to start from
	 * @return the nearest ancestor (or the node itself) of rank {@code SPECIES}, or the original node
	 *         if there is no such ancestor
	 */
	protected static TaxIdNode toSpeciesNode(TaxIdNode node) {
		for (TaxIdNode n = node; n != null; n = n.getParent()) {
			if (Rank.SPECIES.equals(n.getRank())) {
				return n;
			}
		}
		return node;
	}

	/**
	 * Walks up the database's taxonomy from the given node to the nearest ancestor of rank
	 * {@link Rank#SPECIES}. Used to pre-compute the set of expected species keys.
	 *
	 * @param node the node to start from
	 * @return the nearest ancestor (or the node itself) of rank {@code SPECIES}, or the original node
	 *         if there is no such ancestor
	 */
	protected static SmallTaxIdNode toSpeciesNode(SmallTaxIdNode node) {
		for (SmallTaxIdNode n = node; n != null; n = n.getParent()) {
			if (Rank.SPECIES.equals(n.getRank())) {
				return n;
			}
		}
		return node;
	}

	/**
	 * Resolves the well-defined rank to count a node with an indeterminate rank ({@code clade} /
	 * {@code no rank}) under, by walking up to its nearest ranked ancestor.
	 *
	 * @param node the node with an indeterminate rank
	 * @return the rank of the nearest ranked ancestor, or {@code null} if the node has no ranked
	 *         ancestor
	 */
	protected static Rank resolveEffectiveRank(SmallTaxIdNode node) {
		for (SmallTaxIdNode a = node.getParent(); a != null; a = a.getParent()) {
			Rank rank = a.getRank();
			if (rank != null && !rank.isIndeterminate()) {
				return rank;
			}
		}
		return null;
	}

	/**
	 * FASTA reader that, for each k-mer of a region, looks the k-mer up in the database and tallies the
	 * rank of the stored tax id against the region's species in the shared result map.
	 */
	protected class MyFastaReader extends AbstractStoreFastaReader {
		/**
		 * Creates the reader.
		 *
		 * @param bufferSize             the read buffer size
		 * @param taxNodes               the requested tax nodes
		 * @param accessionMap           the accession-to-taxid map
		 * @param k                      the k-mer length
		 * @param maxGenomesPerTaxId     the maximum number of genomes per tax id
		 * @param maxGenomesPerTaxIdRank the rank at which the per-tax-id genome limit applies
		 * @param maxKmersPerTaxId       the maximum number of k-mers per tax id
		 * @param maxDust                the maximum allowed low-complexity (dust) run length
		 * @param stepSize               the k-mer sampling step size
		 * @param completeGenomesOnly    whether to include only complete genomes
		 * @param regionsPerTaxid        the per-taxid region trie
		 * @param enableLowerCaseBases   whether lower-case bases are included
		 */
		public MyFastaReader(int bufferSize, Set<TaxIdNode> taxNodes, AccessionMap accessionMap, int k,
							 int maxGenomesPerTaxId, Rank maxGenomesPerTaxIdRank, long maxKmersPerTaxId, int maxDust,
							 int stepSize, boolean completeGenomesOnly, StringLong2DigitTrie regionsPerTaxid,
							 boolean enableLowerCaseBases) {
			super(bufferSize, taxNodes, accessionMap, k, maxGenomesPerTaxId, maxGenomesPerTaxIdRank, maxKmersPerTaxId,
					maxDust, stepSize, completeGenomesOnly, regionsPerTaxid, enableLowerCaseBases);
		}

		@Override
		protected boolean handleStore() {
			if (node == null || node.getTaxId() == null) {
				return false;
			}
			// Look the k-mer up in the database: which tax id (if any) did it end up under?
			String storedTaxId = store.getLong(byteRingBuffer.getStandardKMer(), null);
			if (storedTaxId == null) {
				// K-mer is not in the database - it did not land anywhere.
				return false;
			}
			SmallTaxIdNode storedNode = taxTree.getNodeByTaxId(storedTaxId);
			if (storedNode == null) {
				return false;
			}
			// Determine the rank slot: a well-defined rank is used directly, an indeterminate one is
			// mapped to the representative ordinal precomputed from the node's tree position.
			int ordinal;
			Rank rank = storedNode.getRank();
			if (rank != null && !rank.isIndeterminate()) {
				ordinal = rank.ordinal();
			} else {
				Integer effectiveOrdinal = effectiveOrdinalByTaxid.get(storedTaxId);
				if (effectiveOrdinal == null) {
					// Indeterminate node with neither a ranked descendant nor a ranked ancestor: it
					// cannot be placed on the rank axis, so it is not counted.
					return false;
				}
				ordinal = effectiveOrdinal;
			}
			// The species node object is shared across reader threads (the tax tree and accession map
			// are shared), so synchronizing on it yields per-species mutual exclusion. Its long[] was
			// pre-created in doMakeThis(), so the map is only read - never structurally modified - here.
			TaxIdNode speciesNode = toSpeciesNode(node);
			long[] counts = map.get(speciesNode.getTaxId());
			if (counts == null) {
				// Should not happen: every tree species was pre-created. Skip rather than modify the
				// map concurrently (which would be unsafe with the fine-grained per-species locking).
				if (getLogger().isWarnEnabled()) {
					getLogger().warn("No pre-created counts for species " + speciesNode.getTaxId() + ".");
				}
				return false;
			}
			synchronized (speciesNode) {
				counts[ordinal]++;
			}
			return true;
		}
	}
}
