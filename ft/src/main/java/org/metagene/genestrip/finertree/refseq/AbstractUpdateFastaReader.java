package org.metagene.genestrip.finertree.refseq;

import org.metagene.genestrip.refseq.AbstractStoreFastaReader;
import org.metagene.genestrip.refseq.AccessionMap;
import org.metagene.genestrip.tax.Rank;
import org.metagene.genestrip.tax.SmallTaxTree;
import org.metagene.genestrip.tax.TaxTree;
import org.metagene.genestrip.util.ByteArrayUtil;

import java.util.Set;

public abstract class AbstractUpdateFastaReader extends AbstractStoreFastaReader {
    private final boolean idNodes;
    private final boolean fileNodes;
    private final boolean dataNodes;

    protected SmallTaxTree.SmallTaxIdNode leafNode;

    public AbstractUpdateFastaReader(int bufferSize, Set<TaxTree.TaxIdNode> taxNodes, AccessionMap accessionMap, int k, int maxGenomesPerTaxId, Rank maxGenomesPerTaxIdRank, long maxKmersPerTaxId, int maxDust, int stepSize, boolean completeGenomesOnly, StringLong2DigitTrie regionsPerTaxid, boolean enableLowerCaseBases, boolean idNodes, boolean fileNodes, boolean dataNodes) {
        super(bufferSize, taxNodes, accessionMap, k, maxGenomesPerTaxId, maxGenomesPerTaxIdRank, maxKmersPerTaxId, maxDust, stepSize, completeGenomesOnly, regionsPerTaxid, enableLowerCaseBases);
        this.idNodes = idNodes;
        this.fileNodes = fileNodes;
        this.dataNodes = dataNodes;
    }

    protected abstract SmallTaxTree getTree();

    @Override
    protected void infoLine() {
        super.infoLine();
        updateLeafNode();
    }

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
