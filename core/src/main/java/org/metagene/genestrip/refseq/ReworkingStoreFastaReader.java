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
package org.metagene.genestrip.refseq;

import java.util.Set;

import org.metagene.genestrip.tax.Rank;
import org.metagene.genestrip.tax.TaxTree;
import org.metagene.genestrip.tax.TaxTree.IDStringGenerator;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;
import org.metagene.genestrip.util.ByteArrayUtil;

/**
 * A store FASTA reader that reworks each contig's tax node into an artificial {@code DATA}/{@code
 * FILE}/{@code ID} child node when the corresponding option is enabled. Shared by the up-front
 * k-mer-counting pass ({@code FillBloomFilterGoal}) and the DB fill ({@code FillDBGoal}) with the same
 * key computation, so both derive identical store values. The counting pass runs in <em>create</em>
 * mode ({@code createNodes = true}): it creates each artificial node in the shared tree if absent
 * (idempotent per {@code (parent, key)} and thread-safe, see {@link TaxTree#dataNode}, so it runs
 * directly in the parallel readers). The fill runs in <em>lookup</em> mode: the counting pass already
 * created every node it needs, so it only looks them up and never mutates the tree.
 */
public abstract class ReworkingStoreFastaReader extends AbstractStoreFastaReader {
	private final TaxTree taxTree;
	private final boolean dataNodes;
	private final boolean fileNodes;
	private final boolean idNodes;
	private final boolean createNodes;
	private final IDStringGenerator idStringGenerator;
	/** Rank a genome is filed at, taxa below it being folded into it; {@code null} to file it where
	 * the taxonomy puts it. See {@code GSConfigKey#FOLD_TAXA_BELOW}. */
	private final Rank foldTaxaBelow;

	/**
	 * Creates the reworking reader.
	 *
	 * @param bufferSize the read buffer size
	 * @param taxNodes the requested tax nodes
	 * @param accessionMap the accession-to-taxid map
	 * @param k the k-mer length
	 * @param maxContigsPerTaxId the maximum number of contigs per tax id
	 * @param maxContigsPerTaxIdRank the rank at which the per-tax-id contig limit applies
	 * @param maxKmersPerTaxId the maximum number of k-mers per tax id
	 * @param maxDust the maximum allowed low-complexity (dust) run length
	 * @param kMerSampling the k-mer sampling step size
	 * @param assemblyAccessionsOnly whether only genomic accessions are considered, dropping `NG_`, `NT_` and `NW_`
	 * @param contigsPerTaxid the per-taxid contig trie
	 * @param enableLowerCaseBases whether lower-case bases are included
	 * @param taxTree the taxonomy tree holding the artificial nodes
	 * @param dataNodes whether to rework into an artificial {@code DATA} node
	 * @param fileNodes whether to rework into an artificial {@code FILE} node
	 * @param idNodes whether to rework into an artificial {@code ID} node
	 * @param createNodes {@code true} to create missing artificial nodes (counting pass), {@code false}
	 *                    to only look up already-created ones (fill)
	 * @param idStringGenerator generator for artificial tax ids (its buffer is mutated, so one per
	 *                          reader); only used (and required) when {@code createNodes} is set
	 * @param foldTaxaBelow the rank a genome is filed at, taxa below it being folded into it, or
	 *                      {@code null} to file it where the taxonomy puts it
	 * @throws IllegalArgumentException if {@code foldTaxaBelow} is set without {@code fileNodes}
	 */
	public ReworkingStoreFastaReader(int bufferSize, Set<TaxIdNode> taxNodes, AccessionMap accessionMap, int k,
			int maxContigsPerTaxId, Rank maxContigsPerTaxIdRank, long maxKmersPerTaxId, int maxDust, int kMerSampling,
			boolean assemblyAccessionsOnly, StringLong2DigitTrie contigsPerTaxid, boolean enableLowerCaseBases,
			TaxTree taxTree, boolean dataNodes, boolean fileNodes, boolean idNodes, boolean createNodes,
			IDStringGenerator idStringGenerator, Rank foldTaxaBelow) {
		super(bufferSize, taxNodes, accessionMap, k, maxContigsPerTaxId, maxContigsPerTaxIdRank, maxKmersPerTaxId,
				maxDust, kMerSampling, assemblyAccessionsOnly, contigsPerTaxid, enableLowerCaseBases);
		checkFoldConfig(foldTaxaBelow, fileNodes);
		this.taxTree = taxTree;
		this.dataNodes = dataNodes;
		this.fileNodes = fileNodes;
		this.idNodes = idNodes;
		this.createNodes = createNodes;
		this.idStringGenerator = idStringGenerator;
		this.foldTaxaBelow = foldTaxaBelow;
	}

	/**
	 * Refuses a fold that is asked for without a file node per genome.
	 * <p>
	 * Checked here and not in the goals, because all three of them build one of these reader and the
	 * combination silently does the opposite of what it is asked for rather than failing: folding
	 * files every genome of a taxon at that taxon, so without a file node of its own each genome
	 * becomes indistinguishable from its neighbours and the taxon is left with nothing below it to
	 * refine. A configuration that reads as {@code file the genomes one rank higher} would then quietly
	 * mean {@code merge them all}, which is worth a refusal rather than a warning.
	 *
	 * @param foldTaxaBelow the rank to file genomes at, or {@code null} for no fold
	 * @param fileNodes whether a file node is created per genome
	 * @throws IllegalArgumentException if a fold is asked for without file nodes
	 */
	static void checkFoldConfig(Rank foldTaxaBelow, boolean fileNodes) {
		if (foldTaxaBelow != null && !fileNodes) {
			throw new IllegalArgumentException("'foldTaxaBelow=" + foldTaxaBelow.getName()
					+ "' requires 'fileNodes=true': folding files every genome of a taxon at the same"
					+ " node, so without a file node of its own each genome becomes indistinguishable"
					+ " from its neighbours and there is nothing left below the taxon to refine.");
		}
	}

	/**
	 * Returns the node a genome resolving to the given one is to be filed at: its lowest ancestor
	 * (itself included) whose rank <em>is</em> {@code below}, or the node unchanged when there is
	 * none and when {@code below} is {@code null}.
	 * <p>
	 * Nodes without a rank are walked through rather than stopped at, and so are nodes whose rank
	 * cannot be ordered against {@code below}: the taxonomy assigns plenty of ranks that the
	 * {@link Rank} enum does not carry -- cohort, parvorder, pathogroup and a dozen more resolve to
	 * {@code null} here -- and stopping at one of them would file a genome at a level nobody asked
	 * for, differently for each lineage that happens to have one.
	 * <p>
	 * The walk stops at the named rank and never passes it. A lineage that carries no node of that
	 * rank -- a rankless node hanging directly under the genus, of which the taxonomy has many --
	 * leaves the genome where it is, and is deliberately not lifted to whatever happens to sit above
	 * instead: asked to file genomes at the species, filing one at its genus would be a coarser
	 * answer than the taxonomy already gave, and one the caller never asked for.
	 * <p>
	 * Nor is a genome ever lifted out of {@code taxNodes}. A contig is only read at all when its
	 * node is one of those (see {@code AbstractRefSeqFastaReader}), so without this check a
	 * {@code taxids.txt} naming a strain, together with a fold at the species, would file that
	 * strain's genomes at a species nobody requested and put a node in the database that is outside
	 * the requested set. Where the target is not requested the genome stays where it is. An empty
	 * {@code taxNodes} means no restriction was asked for and imposes none here either.
	 *
	 * @param node the node the accession resolved to
	 * @param below the rank to file at, or {@code null} to keep the node
	 * @param taxNodes the requested tax nodes, empty for no restriction
	 * @return the node to file the genome at
	 */
	static TaxIdNode foldUp(TaxIdNode node, Rank below, Set<TaxIdNode> taxNodes) {
		if (below == null || node == null) {
			return node;
		}
		for (TaxIdNode n = node; n != null; n = n.getParent()) {
			Rank rank = n.getRank();
			if (rank == null || !rank.isComparableTo(below)) {
				continue;
			}
			if (rank == below) {
				return taxNodes == null || taxNodes.isEmpty() || taxNodes.contains(n) ? n : node;
			}
			if (rank.isAbove(below)) {
				// Past the rank asked for without having met it.
				return node;
			}
		}
		return node;
	}

	@Override
	protected TaxIdNode reworkNode() {
		// The fold happens before anything else and before the node is marked: a taxon nothing is
		// filed at is required by nobody, and markRequired() is the only thing that keeps a node in
		// the SmallTaxTree. The strain nodes a fold skips over therefore disappear by themselves,
		// with no second pass to prune them.
		TaxIdNode res = foldUp(node, foldTaxaBelow, taxNodes);
		res.markRequired();
		if (dataNodes && Rank.DATA.ordinal() != res.getRankOrdinal()) {
			TaxIdNode child = createNodes ? taxTree.dataNode(res, idStringGenerator) : res.getDataChild();
			if (child != null) {
				res = child;
			}
		}
		if (fileNodes && file != null && Rank.FILE.ordinal() != res.getRankOrdinal()) {
			TaxIdNode child = createNodes ? taxTree.fileNode(res, file.getName(), idStringGenerator)
					: res.getChildWithName(file.getName());
			if (child != null) {
				res = child;
			}
		}
		if (idNodes && Rank.ID.ordinal() != res.getRankOrdinal()) {
			int pos = ByteArrayUtil.indexOf(target, 0, size, ' ');
			if (pos < 0) {
				pos = size;
				while (pos > 0 && (target[pos - 1] == '\n' || target[pos - 1] == '\r')) {
					pos--;
				}
			}
			TaxIdNode child = createNodes ? taxTree.idNode(res, target, 1, pos, idStringGenerator)
					: res.getChildWithName(target, 1, pos);
			if (child != null) {
				res = child;
			}
		}
		res.markRequired();
		return res;
	}
}
