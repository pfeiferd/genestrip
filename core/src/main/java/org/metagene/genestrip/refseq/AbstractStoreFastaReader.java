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
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;
import org.metagene.genestrip.util.CGAT;
import org.metagene.genestrip.util.CGATLongBuffer;

/**
 * Abstract RefSeq FASTA reader that streams each region's bases through a {@link CGATLongBuffer}
 * (applying the step size and DUST low-complexity filter) and hands every qualifying k-mer to
 * {@link #handleStore()}.
 */
public abstract class AbstractStoreFastaReader extends AbstractRefSeqFastaReader {
	/** The ring buffer accumulating bases into k-mers. */
	protected final CGATLongBuffer byteRingBuffer;
	/** The number of k-mers dropped as low-complexity (DUST). */
	protected long dustCounter;
	/** The total number of filled k-mers seen. */
	protected long totalKmers;

	private boolean enableLowerCaseBases;

	/**
	 * Creates the reader.
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
	 */
	public AbstractStoreFastaReader(int bufferSize, Set<TaxIdNode> taxNodes, AccessionMap accessionMap, int k, int maxGenomesPerTaxId, Rank maxGenomesPerTaxIdRank,
									long maxKmersPerTaxId, int maxDust, int stepSize, boolean completeGenomesOnly, StringLong2DigitTrie regionsPerTaxid,
									boolean enableLowerCaseBases) {
		super(bufferSize, taxNodes, accessionMap, k, maxGenomesPerTaxId, maxGenomesPerTaxIdRank, maxKmersPerTaxId, stepSize, completeGenomesOnly, regionsPerTaxid);
		byteRingBuffer = new CGATLongBuffer(k, maxDust);
		dustCounter = 0;
		this.enableLowerCaseBases = enableLowerCaseBases;
	}

	/**
	 * Resets the k-mer ring buffer in addition to the superclass region reset.
	 */
	@Override
	protected void startRegion() {
		super.startRegion();
		byteRingBuffer.reset();
	}

	/**
	 * Feeds a data line's bases through the k-mer ring buffer and, every {@code stepSize} bases,
	 * stores each filled non-low-complexity k-mer via {@link #handleStore()}.
	 */
	@Override
	protected void dataLine() {
		if (includeRegion) {
			if (isAllowMoreKmers()) {
				// Strip the trailing line terminator(s) the reader includes in 'size': a single '\n',
				// or '\r\n' for CRLF files. A final line without a trailing newline keeps all bytes
				// (so its last base is not dropped, and a stray '\r' does not reset the ring buffer).
				int end = size;
				while (end > 0 && (target[end - 1] == '\n' || target[end - 1] == '\r')) {
					end--;
				}
				for (int i = 0; i < end; i++) {
					byteRingBuffer.put(enableLowerCaseBases ? CGAT.cgatToUpperCase(target[i]) : target[i]);
					bpsInRegion++;
					if (bpsInRegion % stepSize == 0) {
						if (byteRingBuffer.isFilled()) {
							if (byteRingBuffer.isDust()) {
								dustCounter++;
							} else if (handleStore()) {
								kmersInRegion++;
							}
							totalKmers++;
						}
					}
				}
			}
		}
	}

	/**
	 * Stores the current k-mer (available via the ring buffer); returns whether it was counted as
	 * included.
	 *
	 * @return whether the k-mer was counted as included
	 */
	protected abstract boolean handleStore();
}
