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
 * A store FASTA reader that reworks each region's tax node into an artificial {@code DATA}/{@code
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

	/**
	 * Creates the reworking reader.
	 *
	 * @param bufferSize the read buffer size
	 * @param taxNodes the requested tax nodes
	 * @param accessionMap the accession-to-taxid map
	 * @param k the k-mer length
	 * @param maxGenomesPerTaxId the maximum number of genomes per tax id
	 * @param maxGenomesPerTaxIdRank the rank at which the per-tax-id genome limit applies
	 * @param maxKmersPerTaxId the maximum number of k-mers per tax id
	 * @param maxDust the maximum allowed low-complexity (dust) run length
	 * @param stepSize the k-mer sampling step size
	 * @param completeGenomesOnly whether to include only complete genomes
	 * @param regionsPerTaxid the per-taxid region trie
	 * @param enableLowerCaseBases whether lower-case bases are included
	 * @param taxTree the taxonomy tree holding the artificial nodes
	 * @param dataNodes whether to rework into an artificial {@code DATA} node
	 * @param fileNodes whether to rework into an artificial {@code FILE} node
	 * @param idNodes whether to rework into an artificial {@code ID} node
	 * @param createNodes {@code true} to create missing artificial nodes (counting pass), {@code false}
	 *                    to only look up already-created ones (fill)
	 * @param idStringGenerator generator for artificial tax ids (its buffer is mutated, so one per
	 *                          reader); only used (and required) when {@code createNodes} is set
	 */
	public ReworkingStoreFastaReader(int bufferSize, Set<TaxIdNode> taxNodes, AccessionMap accessionMap, int k,
			int maxGenomesPerTaxId, Rank maxGenomesPerTaxIdRank, long maxKmersPerTaxId, int maxDust, int stepSize,
			boolean completeGenomesOnly, StringLong2DigitTrie regionsPerTaxid, boolean enableLowerCaseBases,
			TaxTree taxTree, boolean dataNodes, boolean fileNodes, boolean idNodes, boolean createNodes,
			IDStringGenerator idStringGenerator) {
		super(bufferSize, taxNodes, accessionMap, k, maxGenomesPerTaxId, maxGenomesPerTaxIdRank, maxKmersPerTaxId,
				maxDust, stepSize, completeGenomesOnly, regionsPerTaxid, enableLowerCaseBases);
		this.taxTree = taxTree;
		this.dataNodes = dataNodes;
		this.fileNodes = fileNodes;
		this.idNodes = idNodes;
		this.createNodes = createNodes;
		this.idStringGenerator = idStringGenerator;
	}

	@Override
	protected TaxIdNode reworkNode() {
		node.markRequired();
		TaxIdNode res = node;
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
