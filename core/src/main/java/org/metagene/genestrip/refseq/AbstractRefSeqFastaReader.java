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

import java.io.File;
import java.io.IOException;
import java.util.Set;

import org.metagene.genestrip.fasta.AbstractFastaReader;
import org.metagene.genestrip.tax.Rank;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;
import org.metagene.genestrip.util.ByteArrayUtil;
import org.metagene.genestrip.util.StringLongDigitTrie;

/**
 * Abstract FASTA reader for RefSeq genome files that resolves each region's tax id via the
 * accession map and tracks, per tax id, the number of regions (genomes) and k-mers already
 * included, so that configurable per-tax-id genome and k-mer limits can be enforced.
 */
public abstract class AbstractRefSeqFastaReader extends AbstractFastaReader {
	/** The set of tax id nodes of interest for this run. */
	protected final Set<TaxIdNode> taxNodes;
	/** Maps sequence accessions to their tax id nodes. */
	protected final AccessionMap accessionMap;
	/** Maximum number of genomes (regions) to include per tax id, or unlimited if non-positive. */
	protected final int maxGenomesPerTaxId;
	/** Maximum number of k-mers to include per tax id, or unlimited if non-positive. */
	protected final long maxKmersPerTaxId;
	/** Taxonomic rank at which the per-tax-id genome limit is applied. */
	protected final Rank maxGenomesPerTaxIdRank;
	/** Per-tax-id trie counting how many regions have already been included. */
	protected final StringLong2DigitTrie regionsPerTaxid;
	/** The k-mer length. */
	protected final int k;
	/** The step size between successive k-mers. */
	protected final int stepSize;
	private final boolean completeGenomesOnly;

	/** Whether the current region is being included. */
	protected boolean includeRegion;
	/** Number of k-mers already included for the current mapped node. */
	protected long kMersForNode;
	/** The tax id node the current region is mapped to. */
	protected TaxIdNode mappedNode;
	/** The tax id node resolved for the current region. */
	protected TaxIdNode node;
	/** The FASTA file currently being read. */
	protected File file;

	/** Whether accession-map lookups are currently bypassed in favor of {@link #mappedNode}. */
	protected boolean ignoreMap;
	/** Number of base pairs seen in the current region. */
	protected long bpsInRegion;
	/** Number of k-mers seen in the current region. */
	protected long kmersInRegion;
	/** Total number of k-mers included so far. */
	protected long includedKmers;


	/**
	 * Creates a RefSeq FASTA reader with the given tax nodes, accession map and per-tax-id limits.
	 *
	 * @param bufferSize             the read buffer size in bytes
	 * @param taxNodes               the set of tax id nodes of interest
	 * @param accessionMap           maps sequence accessions to their tax id nodes
	 * @param k                      the k-mer length
	 * @param maxGenomesPerTaxId     the maximum number of genomes per tax id
	 * @param maxGenomesPerTaxIdRank the rank at which the genome limit is applied
	 * @param maxKmersPerTaxId       the maximum number of k-mers per tax id
	 * @param stepSize               the step size between successive k-mers
	 * @param completeGenomesOnly    whether only complete genomes are considered
	 * @param regionsPerTaxid        the trie counting included regions per tax id
	 */
	public AbstractRefSeqFastaReader(int bufferSize, Set<TaxIdNode> taxNodes, AccessionMap accessionMap, int k, int maxGenomesPerTaxId,
									 Rank maxGenomesPerTaxIdRank, long maxKmersPerTaxId,
									 int stepSize, boolean completeGenomesOnly, StringLong2DigitTrie regionsPerTaxid) {
		super(bufferSize);
		this.taxNodes = taxNodes;
		this.accessionMap = accessionMap;
		this.k = k;
		this.stepSize = stepSize;
		includeRegion = false;
		ignoreMap = false;
		includedKmers = 0;
		this.regionsPerTaxid = regionsPerTaxid;
		this.maxGenomesPerTaxId = maxGenomesPerTaxId;
		this.maxGenomesPerTaxIdRank = maxGenomesPerTaxIdRank;
		this.maxKmersPerTaxId = maxKmersPerTaxId;
		this.completeGenomesOnly = completeGenomesOnly;
	}

	/**
	 * Records the file being read (for reference by subclasses) and delegates to the superclass.
	 */
	public void readFasta(File file) throws IOException {
		this.file = file;
		super.readFasta(file);
	}

	/**
	 * Returns the trie counting how many regions have been included per tax id.
	 *
	 * @return the per-tax-id region-count trie
	 */
	public StringLongDigitTrie getRegionsPerTaxid() {
		return regionsPerTaxid;
	}

	/**
	 * Forces all following regions to be mapped to the given tax id instead of resolving accessions
	 * via the accession map; pass null to re-enable map lookups.
	 *
	 * @param node the tax id node to map all following regions to, or {@code null} to re-enable map lookups
	 */
	public void ignoreAccessionMap(TaxIdNode node) {
		this.ignoreMap = node != null;
		this.mappedNode = node;
	}

	/**
	 * Resets the per-region flags and counters at the start of a new region.
	 */
	@Override
	protected void startRegion() {
		includeRegion = false;
		kmersInRegion = 0;
		bpsInRegion = 0;
	}

	/**
	 * At the end of an included region, adds its k-mer count to its tax id node and all ancestors.
	 */
	@Override
	protected void endRegion() {
		if (includeRegion) {
			includedKmers += kmersInRegion;
			if (node != null) {
				for (TaxIdNode n = node; n != null; n = n.getParent()) {
					regionsPerTaxid.incAndAdd(n.getTaxId(), kmersInRegion);
				}
			}
		}
	}

	/**
	 * Processes a region header: resolves the region's tax id node and decides whether the region is
	 * included, based on the configured per-tax-id genome and k-mer limits.
	 */
	@Override
	protected void infoLine() {
		// Handle new region:
		if (ignoreMap) {
			node = mappedNode;
		}
		else {
			updateNodeFromInfoLine();
		}
		if (node != null && (taxNodes.isEmpty() || taxNodes.contains(node))) {
			node = reworkNode();
			includeRegion = true;
			kMersForNode = 0;
			if (maxGenomesPerTaxIdRank == null) {
				StringLong2DigitTrie.StringLong2 sl = (StringLong2DigitTrie.StringLong2) regionsPerTaxid.get(node.getTaxId());
				if (sl != null) {
					kMersForNode = sl.longValue2;
					if (kMersForNode >= maxKmersPerTaxId || sl.getLongValue() >= maxGenomesPerTaxId) {
						includeRegion = false;
					}
				}
			}
			else {
				for (TaxIdNode n = node; n != null; n = n.getParent()) {
					if (maxGenomesPerTaxIdRank.equals(n.getRank())) {
						StringLong2DigitTrie.StringLong2 sl =
								(StringLong2DigitTrie.StringLong2) regionsPerTaxid.get(n.getTaxId());
						if (sl != null) {
							kMersForNode = sl.longValue2;
							if (kMersForNode >= maxKmersPerTaxId || sl.getLongValue() >= maxGenomesPerTaxId) {
								includeRegion = false;
							}
						}
						break;
					}
				}
			}
		}
		else {
			includeRegion = false;
		}
	}

	/**
	 * Returns whether more k-mers may still be added for the current node without exceeding the
	 * configured per-tax-id k-mer limit.
	 *
	 * @return {@code true} if more k-mers may still be added for the current node
	 */
	public boolean isAllowMoreKmers() {
		return kMersForNode + kmersInRegion < maxKmersPerTaxId;
	}

	/**
	 * Resolves the current region's tax id node from the accession in the header line via the
	 * accession map.
	 */
	protected void updateNodeFromInfoLine() {
		int pos = ByteArrayUtil.indexOf(target, 0, size, ' ');
		if (pos >= 0) {
			node = accessionMap.get(target, 1, pos, completeGenomesOnly);
		}
		else {
			if (getLogger().isWarnEnabled()) {
				getLogger().warn("Inconsistent info line: " + new String(target, 0, size - 1));
			}
			node = null;
		}
	}

	/**
	 * Hook to adjust the resolved tax id node; returns it unchanged by default.
	 *
	 * @return the (possibly adjusted) tax id node for the current region
	 */
	protected TaxIdNode reworkNode() {
		return node;
	}

	/**
	 * A {@link StringLongDigitTrie} whose entries additionally track a second long value (the
	 * accumulated k-mer count) alongside the region count.
	 */
	public static class StringLong2DigitTrie extends StringLongDigitTrie {
		/**
		 * Creates an empty trie.
		 */
		public StringLong2DigitTrie() {
		}

		/**
		 * Increments the region count and adds {@code add} to the k-mer count for the given key,
		 * creating the entry if necessary.
		 *
		 * @param key the tax id key of the entry
		 * @param add the number of k-mers to add to the entry's k-mer count
		 */
		public void incAndAdd(String key, long add) {
			((StringLong2) get(key, this)).incAndAdd(add);
		}

		/**
		 * Creates a new trie entry for the given digit-string key.
		 *
		 * @param digits        the digit-string (tax id) key of the entry
		 * @param createContext the creation context passed through by the trie
		 * @return the new entry
		 */
		@Override
		protected StringLong2 createInGet(String digits, Object createContext) {
			return new StringLong2(digits);
		}

		/**
		 * Creates a new trie entry for the key given as a byte range.
		 *
		 * @param seq           the byte array containing the key
		 * @param start         the start index of the key (inclusive)
		 * @param end           the end index of the key (exclusive)
		 * @param createContext the creation context passed through by the trie
		 * @return the new entry
		 */
		@Override
		protected StringLong createInGet(byte[] seq, int start, int end, Object createContext) {
			return new StringLong2(new String(seq, start, end - start));
		}

		/**
		 * A {@link StringLong} that additionally holds a second long value (the accumulated k-mer
		 * count), where the inherited long value counts regions (genomes).
		 */
		public static class StringLong2 extends StringLong {
			/** The accumulated k-mer count; the inherited long value counts regions (genomes). */
			private long longValue2;

			/**
			 * Creates an entry for the given tax id key with a zero k-mer count.
			 *
			 * @param digits the digit-string (tax id) key of the entry
			 */
			public StringLong2(String digits) {
				super(digits);
			}

			/**
			 * Increments the region count and adds {@code add} to the k-mer count.
			 *
			 * @param add the number of k-mers to add to the k-mer count
			 */
			public synchronized void incAndAdd(long add) {
				longValue++;
				longValue2 += add;
			}

			@Override
			public String toString() {
				return "SL:(taxid:" + stringValue + ", regions: " + longValue + ", kmers: " + longValue2 + ")";
			}
		}
	}
}
