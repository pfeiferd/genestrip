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
import org.metagene.genestrip.util.KMerSampling;

/**
 * Abstract RefSeq FASTA reader that streams each contig's bases through a {@link CGATLongBuffer}
 * (applying the sampling rate and DUST low-complexity filter) and hands every qualifying k-mer to
 * {@link #handleStore(long)}.
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
	 * The threshold selecting one k-mer in {@code kMerSampling}, or {@link KMerSampling#ALL} when every
	 * k-mer is taken. Held as the derived value so that the hot path needs no division.
	 */
	private final long samplingThreshold;

	/**
	 * Creates the reader.
	 *
	 * @param bufferSize the read buffer size
	 * @param taxNodes the requested tax nodes
	 * @param accessionMap the accession-to-taxid map
	 * @param k the k-mer length
	 * @param maxGenomesPerTaxId the maximum number of genomes per tax id
	 * @param maxPerTaxidRank the rank at which the per-tax-id contig limit applies
	 * @param maxKmersPerTaxId the maximum number of k-mers per tax id
	 * @param maxDust the maximum allowed low-complexity (dust) run length
	 * @param kMerSampling the k-mer sampling step size
	 * @param assemblyAccessionsOnly whether only genomic accessions are considered, dropping `NG_`, `NT_` and `NW_`
	 * @param contigsPerTaxid the per-taxid contig trie
	 * @param admittedGenomes the shared set of genome keys admitted so far
	 * @param enableLowerCaseBases whether lower-case bases are included
	 */
	public AbstractStoreFastaReader(int bufferSize, Set<TaxIdNode> taxNodes, AccessionMap accessionMap, int k, int maxGenomesPerTaxId, Rank maxPerTaxidRank,
									long maxKmersPerTaxId, int maxDust, int kMerSampling, boolean assemblyAccessionsOnly, StringLong2DigitTrie contigsPerTaxid, AbstractRefSeqFastaReader.GenomeKeyTrie admittedGenomes,
									boolean enableLowerCaseBases) {
		super(bufferSize, taxNodes, accessionMap, k, maxGenomesPerTaxId, maxPerTaxidRank, maxKmersPerTaxId, kMerSampling, assemblyAccessionsOnly, contigsPerTaxid, admittedGenomes);
		byteRingBuffer = new CGATLongBuffer(k, maxDust);
		dustCounter = 0;
		this.enableLowerCaseBases = enableLowerCaseBases;
		this.samplingThreshold = KMerSampling.thresholdForOneIn(kMerSampling);
	}

	/**
	 * Resets the k-mer ring buffer in addition to the superclass contig reset.
	 */
	@Override
	protected void startContig() {
		super.startContig();
		byteRingBuffer.reset();
	}

	/**
	 * Feeds a data line's bases through the k-mer ring buffer and hands every filled, non-low-complexity
	 * k-mer that belongs to the sample to {@link #handleStore(long)}.
	 */
	@Override
	protected void dataLine() {
		if (includeContig) {
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
					bpsInContig++;
					if (byteRingBuffer.isFilled()) {
						// Which k-mers are kept follows from the k-mer, not from where in the contig it sits,
						// so a k-mer is kept in every genome it occurs in or in none - see KMerSampling for
						// why the database would otherwise misplace the ones it does keep. The canonical
						// encoding travels on to handleStore(), which would only have to compute it again.
						long kmer = byteRingBuffer.getStandardKMer();
						if (KMerSampling.isSampled(kmer, samplingThreshold)) {
							if (byteRingBuffer.isDust()) {
								dustCounter++;
							} else if (handleStore(kmer)) {
								kmersInContig++;
							}
							totalKmers++;
						}
					}
				}
			}
		}
	}

	/**
	 * Stores the given k-mer; returns whether it was counted as included.
	 *
	 * @param kmer the canonical k-mer, as {@link CGATLongBuffer#getStandardKMer()} yields it
	 * @return whether the k-mer was counted as included
	 */
	protected abstract boolean handleStore(long kmer);
}
