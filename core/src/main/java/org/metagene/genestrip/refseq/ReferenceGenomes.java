/*
 * Genestrip
 */
package org.metagene.genestrip.refseq;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.HashSet;
import java.util.Set;

import org.metagene.genestrip.io.BufferedLineReader;

/**
 * The assemblies NCBI marks as a species' reference genome, read off the assembly summary and put
 * in the terms the genome selection works in.
 * <p>
 * The mark is curatorial: it names the assembly to use for a species, and it is what
 * {@code refseq.genomesOnly=preferRefFirst} keeps a place for when a per-taxon limit admits only a
 * few of the genomes on offer. Which assembly carries it has nothing to do with how finished it is
 * -- of the 25,308 marked assemblies of one release only 27.1 per cent are complete -- so a
 * selection that wants the marked one has to recognise drafts as well.
 * <p>
 * Two things are needed for that, and they are found in different places. A complete assembly is
 * already matched to its sequences by {@link AssemblySizeIndex}, so {@link AssemblyInfo#isReference}
 * answers for it. A draft is not in that index, and the join by length would be unsound there
 * anyway; what the summary gives instead is its WGS master accession, whose letter prefix every
 * contig of the assembly carries. RefSeq writes that prefix with {@code NZ_} in front of the
 * GenBank form the summary states -- the summary's {@code AEVD00000000.1} is the catalog's
 * {@code NZ_AEVD00000000.1} and its contigs are {@code NZ_AEVD01000001.1} and so on -- so the
 * prefix, cut the way {@link GenomeKeyTrie} cuts a genome key, names the draft's sequences in the
 * catalog. Checked against one release: of the 89 Streptococcus species whose reference is a draft,
 * 88 are found in the catalog this way, the one miss being a record the release does not carry.
 * <p>
 * This is deliberately not the accession-prefix join that was tried and withdrawn twice for
 * <em>locating assemblies</em>. It does not locate anything: the length-based index stays what
 * finds an assembly, and the prefix is used for the one question the index cannot answer, namely
 * whether a draft is the one its species is represented by.
 */
public class ReferenceGenomes {
	/** Longest summary row taken; the same bound {@link AssemblySizeIndex} uses. */
	private static final int MAX_LINE_SIZE = 8192;

	/** Taxa with a marked reference, by tax id and by species tax id, as the summary states both. */
	private final Set<String> taxaWithReference = new HashSet<String>();
	/** The genome keys of the marked assemblies that are drafts, i.e. that have a WGS master. */
	private final GenomeKeyTrie referenceKeys = new GenomeKeyTrie();
	/** How many of those keys were taken in; the trie counts nothing itself. */
	private int draftCount;

	/**
	 * Reads the summary and collects what the preference needs.
	 *
	 * @param summaryFile the {@code assembly_summary_*.txt} file to read
	 * @throws IOException if the file cannot be read
	 */
	public ReferenceGenomes(File summaryFile) throws IOException {
		byte[] line = new byte[MAX_LINE_SIZE];
		try (BufferedLineReader reader = new BufferedLineReader(new FileInputStream(summaryFile))) {
			int size;
			while ((size = reader.nextLine(line)) > 0) {
				if (size > line.length || line[0] == '#') {
					continue;
				}
				// The four fields wanted, by column: the WGS master, the category, the taxon and the
				// species. Walked once, as the summary is walked everywhere else here.
				int wgsStart = -1, wgsEnd = -1, catStart = -1, catEnd = -1;
				int taxStart = -1, taxEnd = -1, spStart = -1, spEnd = -1;
				int col = 0, from = 0;
				for (int i = 0; i <= size; i++) {
					if (i == size || line[i] == '\t') {
						switch (col) {
						case 3: wgsStart = from; wgsEnd = i; break;
						case 4: catStart = from; catEnd = i; break;
						case 5: taxStart = from; taxEnd = i; break;
						case 6: spStart = from; spEnd = i; break;
						default: break;
						}
						col++;
						from = i + 1;
						if (col > 6) {
							break;
						}
					}
				}
				if (!equals(line, catStart, catEnd, "reference genome")) {
					continue;
				}
				addTaxon(line, taxStart, taxEnd);
				addTaxon(line, spStart, spEnd);
				addKey(line, wgsStart, wgsEnd);
			}
		}
		referenceKeys.freeze();
	}

	/**
	 * Returns whether the summary marks a reference assembly for this taxon, so that a place is
	 * worth keeping for it.
	 *
	 * @param taxId the tax id to ask about, a species where the limit is counted at one
	 * @return whether that taxon has a marked reference assembly
	 */
	public boolean hasReference(String taxId) {
		return taxId != null && taxaWithReference.contains(taxId);
	}

	/**
	 * Returns whether this accession belongs to a draft assembly that is marked as a reference. A
	 * complete one is not recognised here but through {@link AssemblyInfo#isReference()}, which the
	 * caller asks first; the two together cover every marked assembly the catalog carries.
	 *
	 * @param target the buffer holding the accession
	 * @param start the first byte of the accession
	 * @param end the byte after the accession
	 * @return whether the accession's genome is a marked draft reference
	 */
	public boolean isReferenceDraft(byte[] target, int start, int end) {
		return referenceKeys.isAdmitted(target, start, end);
	}

	/**
	 * Returns how many marked assemblies were recognised as drafts, for logging what a build has to
	 * work with.
	 *
	 * @return the number of draft reference genome keys collected
	 */
	public int draftCount() {
		return draftCount;
	}

	/**
	 * Returns the number of taxa a reference is known for, counting a species and the assembly's own
	 * taxon apart, as the summary states them.
	 *
	 * @return the number of tax ids with a marked reference
	 */
	public int taxonCount() {
		return taxaWithReference.size();
	}

	/**
	 * How many genomes a limit node may hold before the genome in hand is turned away: the limit
	 * itself, except that a taxon whose marked assembly is still to come keeps one place back for
	 * it. Without that the limit would be full by the time the marked assembly arrives and which
	 * genomes got in would again follow the order the catalog states them in. The place is held for
	 * the marked assembly alone, so a taxon whose marked assembly the release does not carry ends
	 * one genome short -- the price of holding the limit exactly rather than exceeding it by one.
	 *
	 * @param maxGenomes the configured limit
	 * @param marked whether the genome in hand is the taxon's marked reference
	 * @param referenceAdmitted whether the marked one has already been taken in
	 * @param hasReference whether the summary marks a reference for this taxon at all
	 * @return the number of genomes the taxon may hold for the genome in hand to still fit
	 */
	public static int roomFor(int maxGenomes, boolean marked, boolean referenceAdmitted,
			boolean hasReference) {
		if (marked || referenceAdmitted || !hasReference || maxGenomes == Integer.MAX_VALUE) {
			return maxGenomes;
		}
		return maxGenomes - 1;
	}

	private void addTaxon(byte[] line, int start, int end) {
		if (start >= 0 && end > start) {
			taxaWithReference.add(new String(line, start, end - start));
		}
	}

	/**
	 * Adds the genome key of a marked draft. The summary states the WGS master in its GenBank form
	 * and the catalog carries the RefSeq one, so the key is cut from {@code NZ_} plus that form; a
	 * row without a WGS master is a finished assembly and is left to the length index.
	 */
	private void addKey(byte[] line, int start, int end) {
		if (start < 0 || end <= start || equals(line, start, end, "na")) {
			return;
		}
		byte[] key = new byte[3 + (end - start)];
		key[0] = 'N';
		key[1] = 'Z';
		key[2] = '_';
		System.arraycopy(line, start, key, 3, end - start);
		if (referenceKeys.admit(key, 0, key.length)) {
			draftCount++;
		}
	}

	private static boolean equals(byte[] seq, int start, int end, String word) {
		if (start < 0 || end - start != word.length()) {
			return false;
		}
		for (int i = 0; i < word.length(); i++) {
			if (seq[start + i] != (byte) word.charAt(i)) {
				return false;
			}
		}
		return true;
	}
}
