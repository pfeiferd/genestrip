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
import org.metagene.genestrip.util.DigitTrie;
import org.metagene.genestrip.util.StringLongDigitTrie;

/**
 * Abstract FASTA reader for RefSeq genome files that resolves each contig's tax id via the accession
 * map and tracks, per tax id, the number of contigs and k-mers already included.
 * <p>
 * A contig here is one fasta entry, which is not the same thing as a genome: a draft assembly arrives
 * as one entry per contig, a finished one as one per replicon. Which genomes a database is to hold is
 * therefore not something this reader decides. It asks {@link GenomeKeyTrie}, the selection the
 * accession map carries, and lets a contig in only if its genome is in it - a question of membership,
 * with no counting and no state.
 * <p>
 * That the selection is made elsewhere is the point. Made while reading, it would be made differently
 * in every pass, because which genomes a limit lets in would follow the order the reader threads
 * happen to reach them; a fill that then admitted a genome the counting pass had not seen would file
 * its k-mers at a node the store never registered. Made once from the catalog, it is the same in every
 * pass and in every run. See {@link AccessionMap#getAdmittedGenomes()}.
 * <p>
 * Where no limit is in force the selection is {@code null} and none of this runs.
 */
public abstract class AbstractRefSeqFastaReader extends AbstractFastaReader {
	/** The set of tax id nodes of interest for this run. */
	protected final Set<TaxIdNode> taxNodes;
	/** Maps sequence accessions to their tax id nodes. */
	protected final AccessionMap accessionMap;
	/** Per-tax-id trie counting how many contigs have already been included. */
	protected final StringLong2DigitTrie contigsPerTaxid;
	/**
	 * The genome keys admitted so far, or {@code null} where no genome limit is in force and there is
	 * therefore nothing to remember.
	 * <p>
	 * Handed in like {@link #contigsPerTaxid} and shared with every other reader of the pass - one
	 * reader is created per thread, so an instance of its own would count each genome once per thread.
	 * It must not be static either: a later pass would find every genome of the earlier one admitted
	 * and let them all through.
	 */
	protected final GenomeKeyTrie admittedGenomes;
	/** The k-mer length. */
	protected final int k;
	/** One k-mer in this many is kept, selected by the k-mer itself; see {@code KMerSampling}. */
	protected final int kMerSampling;
	private final boolean assemblyAccessionsOnly;

	/** Whether the current contig is being included. */
	protected boolean includeContig;
	/** The tax id node the current contig is mapped to. */
	protected TaxIdNode mappedNode;
	/** The tax id node resolved for the current contig. */
	protected TaxIdNode node;
	/** The FASTA file currently being read. */
	protected File file;

	/** Whether accession-map lookups are currently bypassed in favor of {@link #mappedNode}. */
	protected boolean ignoreMap;
	/**
	 * Number of base pairs seen in the current contig. Kept as a statistic; it no longer decides which
	 * k-mers are stored, which the k-mer itself does now.
	 */
	protected long bpsInContig;
	/** Number of k-mers seen in the current contig. */
	protected long kmersInContig;
	/** Total number of k-mers included so far. */
	protected long includedKmers;

	/**
	 * Whether the genome limit is in force. Nothing of the genome machinery runs while it is not: no
	 * key is cut from an accession, no trie node is touched and no genome is counted, so a database
	 * that does not ask for the limit pays nothing for it.
	 */
	/** Start of the accession within {@link #target}, or -1 if the info line carried none. */
	private int accessionStart = -1;
	/** End of the accession within {@link #target} (exclusive). */
	private int accessionEnd = -1;


	/**
	 * Creates a RefSeq FASTA reader with the given tax nodes, accession map and per-tax-id limits.
	 *
	 * @param bufferSize             the read buffer size in bytes
	 * @param taxNodes               the set of tax id nodes of interest
	 * @param accessionMap           maps sequence accessions to their tax id nodes
	 * @param k                      the k-mer length
	 * @param kMerSampling               the step size between successive k-mers
	 * @param assemblyAccessionsOnly    whether only genomic accessions are considered, dropping `NG_`, `NT_` and `NW_`
	 * @param contigsPerTaxid        the trie counting included contigs per tax id
	 */
	public AbstractRefSeqFastaReader(int bufferSize, Set<TaxIdNode> taxNodes, AccessionMap accessionMap, int k,
									 int kMerSampling, boolean assemblyAccessionsOnly,
									 StringLong2DigitTrie contigsPerTaxid) {
		super(bufferSize);
		this.taxNodes = taxNodes;
		this.accessionMap = accessionMap;
		this.k = k;
		this.kMerSampling = kMerSampling;
		includeContig = false;
		ignoreMap = false;
		includedKmers = 0;
		this.contigsPerTaxid = contigsPerTaxid;
		// The selection travels with the accession map, which is where it was made - see
		// AccessionMap#getAdmittedGenomes(). Null where no limit is in force, and null where this pass
		// reads no release at all, which is the same thing as far as a cap is concerned.
		this.admittedGenomes = accessionMap == null ? null : accessionMap.getAdmittedGenomes();
		this.assemblyAccessionsOnly = assemblyAccessionsOnly;
	}

	/**
	 * Records the file being read (for reference by subclasses) and delegates to the superclass.
	 */
	public void readFasta(File file) throws IOException {
		this.file = file;
		super.readFasta(file);
	}

	/**
	 * Returns the trie counting how many contigs have been included per tax id.
	 *
	 * @return the per-tax-id contig-count trie
	 */
	public StringLongDigitTrie getContigsPerTaxid() {
		return contigsPerTaxid;
	}

	/**
	 * Forces all following contigs to be mapped to the given tax id instead of resolving accessions
	 * via the accession map; pass null to re-enable map lookups.
	 *
	 * @param node the tax id node to map all following contigs to, or {@code null} to re-enable map lookups
	 */
	public void ignoreAccessionMap(TaxIdNode node) {
		this.ignoreMap = node != null;
		this.mappedNode = node;
	}

	/**
	 * Returns whether the file currently being read is one of the RefSeq release's own fna files.
	 * <p>
	 * The two kinds of input differ in exactly one way, and that is what this reads off: a release
	 * file carries genomes of arbitrary taxa and its contigs are resolved through the accession map,
	 * while a file the project supplies itself - an entry of {@code additional.txt} or a genome
	 * downloaded from Genbank - arrives with the tax node it belongs to and is read with the map
	 * bypassed. {@link #ignoreAccessionMap(TaxIdNode)} is what sets the two apart, and
	 * {@code FastaReaderGoal.readFastas()} calls it before every file so that the answer describes
	 * the file in hand and not one read earlier.
	 * <p>
	 * Callers that treat the release differently from the project's own genomes ask this rather than
	 * {@link #ignoreMap} directly, since it is the distinction and not the map lookup they mean.
	 *
	 * @return whether the current contig stems from the RefSeq release rather than from a project fasta
	 */
	protected final boolean isRefSeqReleaseContig() {
		return !ignoreMap;
	}

	/**
	 * Resets the per-contig flags and counters at the start of a new contig.
	 */
	@Override
	protected void startContig() {
		includeContig = false;
		kmersInContig = 0;
		bpsInContig = 0;
	}

	/**
	 * At the end of an included contig, adds its k-mer count to its tax id node and all ancestors.
	 */
	@Override
	protected void endContig() {
		if (includeContig) {
			includedKmers += kmersInContig;
			if (node != null) {
				for (TaxIdNode n = node; n != null; n = n.getParent()) {
					contigsPerTaxid.incAndAdd(n.getTaxId(), kmersInContig);
				}
			}
		}
	}

	/**
	 * Processes a contig header: resolves the contig's tax id node and decides whether the contig is
	 * included, based on the configured per-tax-id contig and k-mer limits.
	 * <p>
	 * With a {@code maxPerTaxidRank} configured, the limits are counted at the first ancestor of
	 * that rank. A lineage that has no such ancestor - and the taxonomy is full of them - is capped at
	 * its own node instead. Leaving it uncapped, as this once did by falling out of the search loop
	 * without checking anything, means that configuring a rank silently exempts part of the tree from a
	 * limit that was set to bound the database.
	 */
	@Override
	protected void infoLine() {
		// Handle new contig:
		if (ignoreMap) {
			node = mappedNode;
			accessionStart = -1;
			accessionEnd = -1;
		}
		else {
			updateNodeFromInfoLine();
		}
		if (node != null && (taxNodes.isEmpty() || taxNodes.contains(node))) {
			node = reworkNode();
			includeContig = true;
			if (admittedGenomes != null && !isAdmittedGenome()) {
				includeContig = false;
			}
		}
		else {
			includeContig = false;
		}
	}

	/**
	 * Decides whether the genome of the contig in hand is in, counting it at the contig's node and
	 * every ancestor if it is the first contig of it that any reader has seen.
	 * <p>
	 * A genome already admitted stays admitted whatever the counts now say: it entered whole or not at
	 * all, and its remaining contigs are not to be turned away because other genomes filled the limit
	 * meanwhile. That is why the admitted genomes are held rather than the genome in hand - the contigs
	 * of one assembly are scattered through the release files, 204 assemblies arriving in 4,832 runs in
	 * {@code bacteria.1.1.genomic.fna.gz}.
	 * <p>
	 * The limit applies to the RefSeq release only. A fasta the project supplies itself carries no
	 * accession to group by, and it is there because somebody chose it; the release is what needs
	 * taming. Two readers reaching two different genomes at once may both find room for the last place
	 * and take it, which is the same approximation the contig limit makes.
	 *
	 * @param entry the counting entry of the limit node, or null while it has none yet
	 * @return whether the contig's genome is admitted
	 */
	private boolean isAdmittedGenome() {
		// Only the release is capped: a fasta the project supplies itself carries no accession to group
		// by, and it is there because somebody chose it.
		return !isRefSeqReleaseContig() || accessionEnd <= accessionStart
				|| admittedGenomes.isAdmitted(target, accessionStart, accessionEnd);
	}


	/**
	 * Resolves the current contig's tax id node from the accession in the header line via the
	 * accession map.
	 */
	protected void updateNodeFromInfoLine() {
		int pos = ByteArrayUtil.indexOf(target, 0, size, ' ');
		accessionStart = 1;
		accessionEnd = pos;
		if (pos >= 0) {
			node = accessionMap.get(target, 1, pos, assemblyAccessionsOnly);
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
	 * @return the (possibly adjusted) tax id node for the current contig
	 */
	protected TaxIdNode reworkNode() {
		return node;
	}

	/**
	 * A {@link StringLongDigitTrie} whose entries additionally track a second long value (the
	 * accumulated k-mer count) alongside the contig count.
	 */
	public static class StringLong2DigitTrie extends StringLongDigitTrie {
		/**
		 * Creates an empty trie.
		 */
		public StringLong2DigitTrie() {
		}

		/**
		 * Increments the contig count and adds {@code add} to the k-mer count for the given key,
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
		 * A {@link StringLong} that additionally holds the accumulated k-mer count, where the inherited
		 * long value counts contigs.
		 * <p>
		 * Both counters are written under this object's monitor and read without it, by every reader
		 * thread of a pass. A reader deciding whether a contig still fits under a limit may therefore
		 * see counts that are behind, and admit a contig that a stricter accounting would have turned
		 * away: the limits bound the database roughly, not to the assembly. Reading them once all readers
		 * have finished - which is where every other consumer reads them - is safe, since the pass only
		 * ends after their completion has been observed through an atomic counter.
		 */
		public static class StringLong2 extends StringLong {
			/** The accumulated k-mer count; the inherited long value counts contigs. */
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
			 * Increments the contig count and adds {@code add} to the k-mer count.
			 *
			 * @param add the number of k-mers to add to the k-mer count
			 */
			public synchronized void incAndAdd(long add) {
				longValue++;
				longValue2 += add;
			}

			@Override
			public String toString() {
				return "SL:(taxid:" + stringValue + ", contigs: " + longValue + ", kmers: " + longValue2 + ")";
			}
		}
	}
}
