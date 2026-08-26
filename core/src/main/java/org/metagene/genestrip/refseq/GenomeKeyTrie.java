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

import org.metagene.genestrip.util.DigitTrie;

/**
 * The set of genomes a database is to hold, as a trie over the alphabet accessions are written in.
 * <p>
 * A genome is recognised from its accession alone. A WGS accession is a letter prefix naming the
 * sequencing project, a two-digit assembly version and a contig number - {@code NZ_CABEIU010000001} -
 * so every contig of one draft assembly carries the same letter prefix, and
 * {@link #genomeKeyLength(byte[], int, int)} cuts the accession there. A finished replicon has no such
 * prefix and stands for itself, which counts a plasmid, and each segment of a segmented virus, as a
 * genome of its own: an error of one to three per assembly, against the seventy-odd a contig count is
 * out by. Measured against RefSeq release 233 the prefixes recover 8,919 genomes for
 * <em>S. pneumoniae</em> where {@code assembly_summary_refseq.txt} lists 9,263, and 88 to 92 per cent
 * for its neighbours.
 * <p>
 * The set is filled once, while the accession catalog is read, and every pass of a build then reads
 * it. That is what makes {@code maxGenomesPerTaxid} reproducible: deciding while reading would decide
 * differently each time, because which genomes a limit lets in would follow the order the reader
 * threads happen to reach them - and a fill that admitted a genome the counting pass had not seen
 * would file its k-mers at a node the store never registered.
 * <p>
 * A trie because the key is a range of the info line already in hand: looking one up costs a few array
 * reads and allocates nothing, where a hash set would want a {@code String} per contig. Nothing is
 * stored against a key either - the entry's mere presence is the answer, so its value is the cached
 * {@code Boolean.TRUE} and admitting a genome allocates trie nodes and nothing else.
 */
public class GenomeKeyTrie extends DigitTrie<Boolean> {
	private static final long serialVersionUID = 1L;

	/** Passed as the create context of {@link #admit(byte[], int, int)}; the context itself is unused. */
	private static final Object ADMIT = new Object();

	/** Whether the selection is complete and the set is read-only; see {@link #freeze()}. */
	private volatile boolean frozen;

	/**
	 * Creates an empty set of genomes.
	 */
	public GenomeKeyTrie() {
	}

	/**
	 * Returns whether the genome of the accession {@code seq[start, end)} is in the set.
	 *
	 * @param seq   the byte array holding the accession
	 * @param start the start index of the accession (inclusive)
	 * @param end   the end index of the accession (exclusive)
	 * @return whether the accession's genome is admitted
	 */
	public boolean isAdmitted(byte[] seq, int start, int end) {
		return get(seq, start, start + genomeKeyLength(seq, start, end)) != null;
	}

	/**
	 * Adds the genome of the accession {@code seq[start, end)} to the set, and says whether it was not
	 * in it already - which is what makes one genome count once however many accessions it arrives in.
	 *
	 * @param seq   the byte array holding the accession
	 * @param start the start index of the accession (inclusive)
	 * @param end   the end index of the accession (exclusive)
	 * @return whether this call added the genome rather than finding it
	 */
	public boolean admit(byte[] seq, int start, int end) {
		int keyEnd = start + genomeKeyLength(seq, start, end);
		if (get(seq, start, keyEnd) != null) {
			return false;
		}
		get(seq, start, keyEnd, ADMIT);
		return true;
	}

	/**
	 * Closes the set: the selection is complete and nothing may be added to it again.
	 * <p>
	 * Called once, where the catalog scan ends. From there the set is published to every pass of the
	 * build and read by every reader thread without a lock, so a later addition would be both a data
	 * race and a silent change of which genomes the database is to hold - and the passes that already
	 * ran would have read a different selection from the ones that follow. That is the failure this
	 * design exists to prevent, so it is refused rather than trusted not to happen.
	 */
	public void freeze() {
		frozen = true;
	}

	/**
	 * Returns whether the selection is closed.
	 *
	 * @return whether nothing may be added any more
	 */
	public boolean isFrozen() {
		return frozen;
	}

	/** Digits, then the upper-case letters, then the underscore, then anything else. */
	@Override
	protected int mapToIndex(byte bite, int pos) {
		if (bite >= '0' && bite <= '9') {
			return bite - '0';
		}
		if (bite >= 'A' && bite <= 'Z') {
			return bite - 'A' + 10;
		}
		return bite == '_' ? 36 : 37;
	}

	@Override
	protected int range(int pos) {
		return 38;
	}

	/**
	 * Marks the genome as admitted - the entry's presence is the whole of it, so the value is the
	 * cached {@code TRUE} and nothing is allocated for it.
	 * <p>
	 * This is also where a closed selection refuses an addition. It is checked here rather than in
	 * {@link #admit(byte[], int, int)} because this is the one place a key can come into being, so a
	 * caller reaching past that method through {@code DigitTrie} is refused too.
	 *
	 * @param seq           the byte array holding the genome key
	 * @param start         the start index of the key (inclusive)
	 * @param end           the end index of the key (exclusive)
	 * @param createContext ignored
	 * @return {@code TRUE}
	 */
	@Override
	protected Boolean createInGet(byte[] seq, int start, int end, Object createContext) {
		if (frozen) {
			throw new IllegalStateException("The genome selection is complete and cannot take "
					+ new String(seq, start, end - start) + ": it is read by every pass of the build, and"
					+ " adding to it would change which genomes the database holds halfway through.");
		}
		return Boolean.TRUE;
	}

	/**
	 * Returns how much of the given accession names the genome it belongs to.
	 * <p>
	 * The assembly version is deliberately left out of the key: {@code CABEIU01} and {@code CABEIU02}
	 * are two versions of one assembly, and counting them apart doubles the tally - against RefSeq
	 * release 233 it turns 8,919 genomes of <em>S. pneumoniae</em> into 16,771. Anything that is not a
	 * WGS accession names itself, minus its version.
	 *
	 * @param seq   the byte array holding the accession
	 * @param start the start index of the accession (inclusive)
	 * @param end   the end index of the accession (exclusive)
	 * @return the length of the leading part of the accession that identifies its genome
	 */
	public static int genomeKeyLength(byte[] seq, int start, int end) {
		int stop = start;
		while (stop < end && seq[stop] != '.') {
			stop++;
		}
		int i = start;
		while (i < stop && seq[i] != '_') {
			i++;
		}
		if (i == stop) {
			return stop - start;
		}
		int lettersStart = ++i;
		while (i < stop && seq[i] >= 'A' && seq[i] <= 'Z') {
			i++;
		}
		int digits = i;
		while (digits < stop && seq[digits] >= '0' && seq[digits] <= '9') {
			digits++;
		}
		if (i - lettersStart >= 4 && digits - i >= 8 && digits == stop) {
			return i - start;
		}
		return stop - start;
	}
}
