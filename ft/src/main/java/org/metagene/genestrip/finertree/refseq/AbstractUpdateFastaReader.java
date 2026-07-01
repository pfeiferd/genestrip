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
package org.metagene.genestrip.finertree.refseq;

import org.metagene.genestrip.refseq.AbstractStoreFastaReader;
import org.metagene.genestrip.refseq.AccessionMap;
import org.metagene.genestrip.tax.Rank;
import org.metagene.genestrip.tax.SmallTaxTree;
import org.metagene.genestrip.tax.TaxTree;
import org.metagene.genestrip.util.ByteArrayUtil;

import java.util.Set;

/**
 * Base fasta reader that, for each info line, resolves the leaf taxonomy node to which the following
 * *k*-mers belong, optionally descending into artificial id, file or data nodes of the finer tree.
 */
public abstract class AbstractUpdateFastaReader extends AbstractStoreFastaReader {
    private final boolean idNodes;
    private final boolean fileNodes;
    private final boolean dataNodes;

    /** The leaf node of the finer tree to which the k-mers of the current region belong. */
    protected SmallTaxTree.SmallTaxIdNode leafNode;

    /**
     * Creates a reader that resolves finer-tree leaf nodes for the k-mers of each region.
     *
     * @param bufferSize             the input read buffer size
     * @param taxNodes               the tax nodes to be included
     * @param accessionMap           the map from accession numbers to tax nodes
     * @param k                      the k-mer length
     * @param maxGenomesPerTaxId     the maximum number of genomes to consider per tax id
     * @param maxGenomesPerTaxIdRank the rank at which the per-tax-id genome limit applies
     * @param maxKmersPerTaxId       the maximum number of k-mers to store per tax id
     * @param maxDust                the maximum dust (low-complexity) threshold
     * @param stepSize               the step size between stored k-mers
     * @param completeGenomesOnly    whether only complete genomes are considered
     * @param regionsPerTaxid        the trie counting regions per tax id
     * @param enableLowerCaseBases   whether lower-case bases are treated as regular bases
     * @param idNodes                whether to descend into artificial id nodes
     * @param fileNodes              whether to descend into artificial file nodes
     * @param dataNodes              whether to descend into artificial data nodes
     */
    public AbstractUpdateFastaReader(int bufferSize, Set<TaxTree.TaxIdNode> taxNodes, AccessionMap accessionMap, int k, int maxGenomesPerTaxId, Rank maxGenomesPerTaxIdRank, long maxKmersPerTaxId, int maxDust, int stepSize, boolean completeGenomesOnly, StringLong2DigitTrie regionsPerTaxid, boolean enableLowerCaseBases, boolean idNodes, boolean fileNodes, boolean dataNodes) {
        super(bufferSize, taxNodes, accessionMap, k, maxGenomesPerTaxId, maxGenomesPerTaxIdRank, maxKmersPerTaxId, maxDust, stepSize, completeGenomesOnly, regionsPerTaxid, enableLowerCaseBases);
        this.idNodes = idNodes;
        this.fileNodes = fileNodes;
        this.dataNodes = dataNodes;
    }

    /**
     * Returns the taxonomy tree used to resolve finer-tree nodes.
     *
     * @return the taxonomy tree used to resolve nodes for read tax ids
     */
    protected abstract SmallTaxTree getTree();

    @Override
    protected void infoLine() {
        super.infoLine();
        updateLeafNode();
    }

    /**
     * Sets {@link #leafNode} to the node for the current region's tax id, descending into a matching
     * id, file or data child node when the respective feature is enabled.
     */
    protected void updateLeafNode() {
        if (includeRegion) {
            leafNode = getTree().getNodeByTaxId(node.getTaxId());
            if (leafNode != null) {
                // Try to find id node first:
                // If we first try the file node via the file's name we may end up at the wrong node...
                if (idNodes) {
                    if (Rank.ID.ordinal() != leafNode.getRankOrdinal()) {
                        int pos = ByteArrayUtil.indexOf(target, 0, size, ' ');
                        if (pos < 0) {
                            pos = size;
                        }
                        SmallTaxTree.SmallTaxIdNode h = leafNode.getDescendantWithName(target, 1, pos);
                        if (h != null) {
                            leafNode = h;
                            return;
                        }
                    }
                }
                if (fileNodes && file != null) {
                    if (Rank.FILE.ordinal() != leafNode.getRankOrdinal()) {
                        SmallTaxTree.SmallTaxIdNode h = leafNode.getDescendantWithName(file.getName());
                        if (h != null) {
                            leafNode = h;
                            return;
                        }
                    }
                }
                if (dataNodes) {
                    if (Rank.DATA.ordinal() != leafNode.getRankOrdinal()) {
                        SmallTaxTree.SmallTaxIdNode h = leafNode.getDataChild();
                        if (h != null) {
                            leafNode = h;
                            return;
                        }
                    }
                }
            }
        }
    }
}
