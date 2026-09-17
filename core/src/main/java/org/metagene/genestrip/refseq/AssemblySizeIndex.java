/*
 * Genestrip
 */
package org.metagene.genestrip.refseq;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.EnumSet;
import java.util.HashSet;
import java.util.Set;

import org.metagene.genestrip.io.BufferedLineReader;

import org.metagene.genestrip.genbank.AssemblySummaryReader.AssemblyQuality;

/**
 * An index over NCBI's assembly summary, keyed by the taxon and the assembly's total length.
 * <p>
 * The accession catalog a database is built from names sequences, not assemblies, and carries no
 * assembly accession, so the summary cannot be joined to it directly. What the two do share is a
 * taxon and a number of bases: an assembly's {@code genome_size} is the sum of the lengths of the
 * sequences it consists of. Sorting the summary by that pair and searching it for a sum taken from
 * the catalog therefore locates the assembly.
 * <p>
 * Only assemblies whose length can be taken for the organism's are indexed - complete and
 * chromosome-level ones, {@link #COMPLETE_OR_CHROMOSOME} - which is the question the index is there
 * to answer: a hit says the sum is that of such an assembly and a miss says it is not. Screening the
 * level after the search instead would lose a sum whose admitted assembly lay within the gap but
 * which a draft happened to lie nearer to. Four assemblies in five being drafts, leaving them out
 * also makes the index a fraction of the size.
 * <p>
 * The match is not required to be exact. A submitted total and the sum of the released sequences
 * differ by a few bases often enough to matter, so a search admits any entry within a gap of the
 * sum and takes the nearest. Measured over <em>E. coli</em>,
 * <em>K. pneumoniae</em>, <em>V. cholerae</em> and <em>P. aeruginosa</em>, an exact search locates
 * 74.8 per cent of the complete assemblies and a gap of ten locates 82.1 per cent; a gap of a
 * hundred reaches 92.0 per cent but begins to group sequences that do not belong together, which
 * shows as more assemblies being 'located' than exist.
 * <p>
 * Both size columns are indexed, {@code genome_size} and {@code genome_size_ungapped}, because which
 * of the two a sum of released sequences equals varies: the ungapped column matches 96.7 per cent of
 * the whole-genome shotgun master records against 65.0 per cent for the gapped one.
 */
public class AssemblySizeIndex {

	/** Longest summary row taken; rows are a few hundred bytes and the longest seen is far below this. */
	private static final int MAX_LINE_SIZE = 8192;

	/**
	 * The taxon of each entry, sorted by taxon and then by size. Held beside {@link #sizes} rather
	 * than packed with it into one long: indexing complete assemblies only leaves a few hundred
	 * thousand entries, so the four bytes a plain taxon costs over a packed one are a few megabytes
	 * altogether, and a size does not have to be squeezed into the bits a taxon leaves over. The
	 * largest genome in the summary is forty gigabases and would not fit in thirty-two of them.
	 */
	private final int[] taxids;

	/** The total length of each entry's assembly, in the same order. */
	private final long[] sizes;

	/** Whether the summary marks each entry's assembly as its species' reference genome. */
	private final boolean[] references;

	/**
	 * Which assembly each entry belongs to, in the same order.
	 * <p>
	 * An assembly is indexed under both its taxon and its species and under both size columns, so it
	 * owns up to four entries; they share one id. It stands for the assembly accession that the
	 * accession catalog does not carry, and lets a caller tell whether two sums it matched are the
	 * same genome found twice - which they often are, because a genome's chromosome alone can match
	 * a different assembly of the species more closely than the whole genome matches its own.
	 */
	private final int[] assemblyIds;

	/** How many assemblies are indexed, which is one more than the highest id. */
	private final int assemblyCount;

	/**
	 * The assembly level of each entry, in the same order, held as the enum's ordinal in one byte.
	 * Every entry is of an admitted level, so this only distinguishes within those.
	 */
	private final byte[] levels;

	/** The enum constants, cached so that a lookup does not clone the array each time. */
	private static final AssemblyQuality[] QUALITIES = AssemblyQuality.values();

	/**
	 * The levels indexed unless a caller says otherwise, and the same ones the Genbank path admits
	 * by default through {@code genbank.fastaQualities}.
	 * <p>
	 * A complete assembly is gapless by definition - of the summary's seventy-three thousand of
	 * them, not one has a {@code genome_size} differing from its {@code genome_size_ungapped}. A
	 * chromosome-level assembly falls about one per cent short of the complete genomes of its own
	 * species, which is little enough to admit it. That one per cent is not what the two size columns
	 * say: their ratio counts only the gaps the assembler declared, a median of 0.02 per cent against
	 * the 1.18 per cent the same-species comparison implies, so the columns are not a completeness and
	 * are not used as one here. A scaffold or contig assembly does not represent what it is missing at
	 * all, which is why those are left out and why a contig assembly's gap fraction reads as zero.
	 */
	public static final Set<AssemblyQuality> COMPLETE_OR_CHROMOSOME = Collections.unmodifiableSet(EnumSet.of(
			AssemblyQuality.COMPLETE_LATEST, AssemblyQuality.COMPLETE, AssemblyQuality.CHROMOSOME_LATEST,
			AssemblyQuality.CHROMOSOME));

	/** Complete assemblies only, which are gapless by definition. */
	public static final Set<AssemblyQuality> COMPLETE_ONLY = Collections
			.unmodifiableSet(EnumSet.of(AssemblyQuality.COMPLETE_LATEST, AssemblyQuality.COMPLETE));

	/** Chromosome-level assemblies, which usually carry gaps of about one per cent of the genome. */
	public static final Set<AssemblyQuality> CHROMOSOME_LEVEL = Collections
			.unmodifiableSet(EnumSet.of(AssemblyQuality.CHROMOSOME_LATEST, AssemblyQuality.CHROMOSOME));

	/** Every level, for a caller wanting the summary's metadata rather than an answer about completeness. */
	public static final Set<AssemblyQuality> ALL_LEVELS = Collections
			.unmodifiableSet(EnumSet.allOf(AssemblyQuality.class));

	/**
	 * Reads the summary file and builds the index.
	 *
	 * Only {@link #COMPLETE_OR_CHROMOSOME} is indexed.
	 *
	 * @param summaryFile the {@code assembly_summary_*.txt} file to read
	 * @throws IOException if the file cannot be read
	 */
	public AssemblySizeIndex(File summaryFile) throws IOException {
		this(summaryFile, COMPLETE_OR_CHROMOSOME);
	}

	/**
	 * Reads the summary file and builds the index over the given levels.
	 * <p>
	 * Which levels are indexed is which levels a search can return: a hit says the sum is that of an
	 * admitted assembly and a miss says it is not. Screening the level after the search instead
	 * would lose a sum whose admitted assembly lay within the gap but which an unadmitted one
	 * happened to lie nearer to - at the default gap of a hundred, <em>S. aureus</em> already has a
	 * contig assembly three bases from a sum whose complete assembly is seventy-eight away.
	 *
	 * @param summaryFile the {@code assembly_summary_*.txt} file to read
	 * @param admitted the levels to index; {@link #ALL_LEVELS} for all of them, in which case
	 *            {@link Match#getLevel()} is what says which level was found
	 * @throws IOException if the file cannot be read
	 */
	public AssemblySizeIndex(File summaryFile, Collection<AssemblyQuality> admitted) throws IOException {
		this(summaryFile, admitted, false);
	}

	/**
	 * Reads the summary file and builds the index over the given levels, optionally keeping only the
	 * assembly NCBI marks as a species' reference where it has marked one.
	 *
	 * @param summaryFile the {@code assembly_summary_*.txt} file to read
	 * @param admitted the levels to index
	 * @param preferReference whether a species that has an admitted assembly marked
	 *            {@code reference genome} should be represented by that assembly alone. A species
	 *            with no such mark keeps all of its admitted assemblies. This chooses among
	 *            assemblies and does not find them: most reference genomes are not complete, so the
	 *            mark is useless as a level filter and is only good for breaking a tie the level
	 *            filter has already left open.
	 * @throws IOException if the file cannot be read
	 */
	public AssemblySizeIndex(File summaryFile, Collection<AssemblyQuality> admitted, boolean preferReference)
			throws IOException {
		// Tested per row, so the set is flattened to a lookup by ordinal.
		boolean[] admit = new boolean[QUALITIES.length];
		for (AssemblyQuality q : admitted) {
			admit[q.ordinal()] = true;
		}
		int[] t = new int[1 << 16];
		long[] z = new long[1 << 16];
		int[] ids = new int[1 << 16];
		int assembly = 0;
		// Per admitted assembly, what the reference preference needs: the species it belongs to and
		// whether the summary marks it. Grown alongside the entry arrays but indexed by assembly.
		int[] aSpecies = new int[1 << 12];
		boolean[] aRef = new boolean[1 << 12];
		byte[] q = new byte[1 << 16];
		int n = 0;
		// Parsed off the bytes rather than through split(): the file has half a million rows of
		// thirty-eight columns, so splitting would make twenty million strings to read five fields
		// from, and this is the same reason the accession catalog is read the way it is.
		byte[] line = new byte[MAX_LINE_SIZE];
		try (BufferedLineReader reader = new BufferedLineReader(new FileInputStream(summaryFile))) {
			int size;
			while ((size = reader.nextLine(line)) > 0) {
				if (size > line.length || line[0] == '#') {
					continue;
				}
				// The five fields wanted, by column: the taxon, the species, the level, and the two
				// size columns. Walking the tabs once is cheaper than indexing to each separately.
				int col = 0;
				int catStart = -1, catEnd = -1;
				int taxStart = -1, taxEnd = -1, spStart = -1, spEnd = -1;
				int lvlStart = -1, lvlEnd = -1, verStart = -1, verEnd = -1;
				int gsStart = -1, gsEnd = -1, guStart = -1, guEnd = -1;
				int from = 0;
				for (int i = 0; i <= size; i++) {
					if (i == size || line[i] == '\t') {
						switch (col) {
						case 4: catStart = from; catEnd = i; break;
						case 5: taxStart = from; taxEnd = i; break;
						case 6: spStart = from; spEnd = i; break;
						case 10: verStart = from; verEnd = i; break;
						case 11: lvlStart = from; lvlEnd = i; break;
						case 25: gsStart = from; gsEnd = i; break;
						case 26: guStart = from; guEnd = i; break;
						default: break;
						}
						col++;
						from = i + 1;
						if (col > 26) {
							break;
						}
					}
				}
				if (gsStart < 0 && guStart < 0) {
					continue;
				}
				byte ord = (byte) levelOrdinal(line, lvlStart, lvlEnd, verStart, verEnd);
				if (!admit[ord]) {
					continue;
				}
				// One id per admitted row, shared by the up-to-four entries the row is indexed under.
				int id = assembly++;
				if (id == aSpecies.length) {
					aSpecies = Arrays.copyOf(aSpecies, aSpecies.length * 2);
					aRef = Arrays.copyOf(aRef, aRef.length * 2);
				}
				aSpecies[id] = (int) parseLong(line, spStart, spEnd);
				aRef[id] = equals(line, catStart, catEnd, "reference genome");
				// Both the taxon and the species are indexed: the catalog files a sequence under
				// whichever of the two the submitter used, and it is not always the one the summary
				// names first.
				for (int ti = 0; ti < 2; ti++) {
					long tax = ti == 0 ? parseLong(line, taxStart, taxEnd) : parseLong(line, spStart, spEnd);
					if (tax < 0) {
						continue;
					}
					for (int si = 0; si < 2; si++) {
						long sz = si == 0 ? parseLong(line, gsStart, gsEnd) : parseLong(line, guStart, guEnd);
						if (sz <= 0) {
							continue;
						}
						if (n == t.length) {
							t = Arrays.copyOf(t, t.length * 2);
							z = Arrays.copyOf(z, z.length * 2);
							ids = Arrays.copyOf(ids, ids.length * 2);
							q = Arrays.copyOf(q, q.length * 2);
						}
						t[n] = (int) tax;
						z[n] = sz;
						ids[n] = id;
						q[n] = ord;
						n++;
					}
				}
			}
		}
		// Where a reference is preferred, a species that has one marked is represented by that
		// assembly alone and its other assemblies are dropped here, before anything is laid out. The
		// mark is on the assembly and the preference is over the species, so the species that have
		// one have to be collected before any assembly can be judged.
		if (preferReference) {
			Set<Integer> speciesWithReference = new HashSet<Integer>();
			for (int i = 0; i < assembly; i++) {
				if (aRef[i]) {
					speciesWithReference.add(aSpecies[i]);
				}
			}
			int[] remap = new int[assembly];
			int kept = 0;
			for (int i = 0; i < assembly; i++) {
				remap[i] = aRef[i] || !speciesWithReference.contains(aSpecies[i]) ? kept++ : -1;
			}
			int m = 0;
			for (int i = 0; i < n; i++) {
				int to = remap[ids[i]];
				if (to >= 0) {
					t[m] = t[i];
					z[m] = z[i];
					ids[m] = to;
					q[m] = q[i];
					m++;
				}
			}
			boolean[] rr = new boolean[kept];
			for (int i = 0; i < assembly; i++) {
				if (remap[i] >= 0) {
					rr[remap[i]] = aRef[i];
				}
			}
			aRef = rr;
			n = m;
			assembly = kept;
		}

		// Sorting the parallel arrays: Arrays.sort cannot carry them along, so an order is sorted
		// instead and the arrays are laid out afterwards.
		Integer[] order = new Integer[n];
		for (int i = 0; i < n; i++) {
			order[i] = i;
		}
		final int[] tt = t;
		final long[] zz = z;
		Arrays.sort(order, (a, b) -> tt[a] != tt[b] ? Integer.compare(tt[a], tt[b]) : Long.compare(zz[a], zz[b]));
		taxids = new int[n];
		sizes = new long[n];
		assemblyIds = new int[n];
		levels = new byte[n];
		references = new boolean[n];
		for (int i = 0; i < n; i++) {
			taxids[i] = t[order[i]];
			sizes[i] = z[order[i]];
			assemblyIds[i] = ids[order[i]];
			levels[i] = q[order[i]];
			references[i] = aRef[assemblyIds[i]];
		}
		assemblyCount = assembly;
	}

	/**
	 * Returns how many assemblies are indexed, which is what an id is bounded by.
	 *
	 * @return the number of assemblies indexed
	 */
	public int getAssemblyCount() {
		return assemblyCount;
	}

	/**
	 * Returns the ordinal of the assembly level of a row, read off the bytes.
	 * <p>
	 * {@link AssemblyQuality#fromString} would do this given two strings, and two strings per row is a
	 * million of them over a summary, for five fields that are each one of a handful of fixed words.
	 * The words are compared where they lie instead. The mapping is the one that method makes, and has
	 * to be kept beside it.
	 */
	private static int levelOrdinal(byte[] seq, int lvlStart, int lvlEnd, int verStart, int verEnd) {
		boolean latest = equals(seq, verStart, verEnd, "latest");
		if (equals(seq, lvlStart, lvlEnd, "Complete Genome")) {
			return (latest ? AssemblyQuality.COMPLETE_LATEST : AssemblyQuality.COMPLETE).ordinal();
		}
		if (equals(seq, lvlStart, lvlEnd, "Chromosome")) {
			return (latest ? AssemblyQuality.CHROMOSOME_LATEST : AssemblyQuality.CHROMOSOME).ordinal();
		}
		if (equals(seq, lvlStart, lvlEnd, "Scaffold")) {
			return (latest ? AssemblyQuality.SCAFFOLD_LATEST : AssemblyQuality.SCAFFOLD).ordinal();
		}
		if (equals(seq, lvlStart, lvlEnd, "Contig")) {
			return (latest ? AssemblyQuality.CONTIG_LATEST : AssemblyQuality.CONTIG).ordinal();
		}
		return (latest ? AssemblyQuality.LATEST : AssemblyQuality.NONE).ordinal();
	}

	/**
	 * Returns whether {@code seq[start, end)} is the given word. The word is a constant of this class,
	 * so nothing is allocated to make the comparison.
	 */
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

	/** Parses a non-negative long from bytes, returning -1 where the field is not one. */
	private static long parseLong(byte[] seq, int start, int end) {
		if (start < 0 || start >= end) {
			return -1;
		}
		long v = 0;
		for (int i = start; i < end; i++) {
			byte c = seq[i];
			if (c < '0' || c > '9') {
				return -1;
			}
			v = v * 10 + (c - '0');
		}
		return v;
	}

	/** What a search found: the level of the nearest assembly and how far off its length was. */
	public static class Match {
		private final AssemblyQuality level;
		private final long distance;
		private final int assemblyId;
		private final boolean reference;

		Match(AssemblyQuality level, long distance, int assemblyId, boolean reference) {
			this.level = level;
			this.distance = distance;
			this.assemblyId = assemblyId;
			this.reference = reference;
		}

		/**
		 * Returns whether the summary marks this assembly as its species' reference genome, i.e. as
		 * NCBI's own choice of the assembly to use for the species. It says nothing about how
		 * finished the assembly is: most reference genomes are not complete.
		 *
		 * @return whether the assembly is marked as a reference genome
		 */
		public boolean isReference() {
			return reference;
		}

		/**
		 * Returns which assembly was matched, as an id below
		 * {@link AssemblySizeIndex#getAssemblyCount()}.
		 * <p>
		 * It stands in for the assembly accession, which the accession catalog does not carry. Two
		 * matches with the same id are the same genome, so a caller counting genomes can tell a
		 * second sighting from a second genome.
		 *
		 * @return the id of the assembly matched
		 */
		public int getAssemblyId() {
			return assemblyId;
		}

		/**
		 * Returns the assembly level of the nearest assembly within the gap.
		 *
		 * @return the assembly level of the nearest assembly within the gap
		 */
		public AssemblyQuality getLevel() {
			return level;
		}

		/**
		 * Returns how far the sum was from the assembly's recorded length.
		 *
		 * @return how many bases the sum was from that assembly's recorded length
		 */
		public long getDistance() {
			return distance;
		}

		/**
		 * Returns whether the assembly is complete, and so gapless: no complete assembly in the
		 * summary has a {@code genome_size} differing from its {@code genome_size_ungapped}.
		 *
		 * @return whether the assembly is complete
		 */
		public boolean isComplete() {
			return COMPLETE_ONLY.contains(level);
		}

		/**
		 * Returns whether the assembly is chromosome-level rather than complete. Nine in ten such
		 * assemblies carry gaps, averaging about one per cent of the genome, so a caller reckoning
		 * with how much of an organism the database holds wants to know which of the two it has.
		 *
		 * @return whether the assembly is chromosome-level
		 */
		public boolean isChromosome() {
			return CHROMOSOME_LEVEL.contains(level);
		}
	}

	/**
	 * Returns the nearest assembly of the given taxon whose recorded length is within {@code gap} of
	 * the given sum, or null where there is none. Where several are within the gap the nearest wins,
	 * which is what lets a caller comparing several groupings of the same sequences prefer the one
	 * that fits best.
	 *
	 * @param taxid the taxon the sequences were filed under
	 * @param size the sum of the lengths of the sequences taken to be one assembly
	 * @param gap how far the sum may be from an assembly's recorded length
	 * @return the nearest match, or null
	 */
	public Match nearest(int taxid, long size, int gap) {
		// The entries of a taxon are sorted by size, so the nearest to the sum is either the first
		// entry at or above it or the last one below it - its two neighbours. Searching for the sum
		// itself and looking at both is therefore the whole of the work, and the gap only decides
		// afterwards whether the nearer of the two is near enough. Scanning the window the gap opens
		// would find the same entry and take as long as the window is wide.
		int i = firstAtLeast(taxid, size);
		int best = -1;
		long bestDist = Long.MAX_VALUE;
		int bestId = -1;
		boolean bestRef = false;
		if (i < taxids.length && taxids[i] == taxid && sizes[i] - size <= gap) {
			bestDist = sizes[i] - size;
			best = levels[i];
			bestId = assemblyIds[i];
			bestRef = references[i];
		}
		if (i > 0 && taxids[i - 1] == taxid) {
			long d = size - sizes[i - 1];
			// Strictly nearer, so a tie goes to the entry at or above the sum. Which of two equally
			// distant assemblies is returned does not matter: the caller asks how far off the sum
			// is, and both answer the same.
			if (d <= gap && d < bestDist) {
				bestDist = d;
				best = levels[i - 1];
				bestId = assemblyIds[i - 1];
				bestRef = references[i - 1];
			}
		}
		return best < 0 ? null : new Match(QUALITIES[best], bestDist, bestId, bestRef);
	}

	/**
	 * Returns the distance of the nearest entry found by scanning the whole window the gap opens,
	 * or null where the window holds none. This is what {@link #nearest(int, long, int)} did before
	 * it was reduced to looking at the two neighbours of the sum, and it is kept as the reference
	 * the test compares that against: the reduction rests on the entries of a taxon being sorted by
	 * size, and an argument of that kind is worth checking against the thing it replaced.
	 *
	 * @param taxid the taxon the sequences were filed under
	 * @param size the sum of the lengths of the sequences taken to be one assembly
	 * @param gap how far the sum may be from an assembly's recorded length
	 * @return the distance of the nearest entry, or null
	 */
	Long nearestByScan(int taxid, long size, int gap) {
		long bestDist = Long.MAX_VALUE;
		for (int j = firstAtLeast(taxid, Math.max(0, size - gap)); j < taxids.length; j++) {
			if (taxids[j] != taxid || sizes[j] > size + gap) {
				break;
			}
			bestDist = Math.min(bestDist, Math.abs(sizes[j] - size));
		}
		return bestDist == Long.MAX_VALUE ? null : bestDist;
	}

	/** Returns the first index whose taxon and size are not below the given pair. */
	private int firstAtLeast(int taxid, long size) {
		int lo = 0;
		int hi = taxids.length;
		while (lo < hi) {
			int mid = (lo + hi) >>> 1;
			if (taxids[mid] < taxid || (taxids[mid] == taxid && sizes[mid] < size)) {
				lo = mid + 1;
			} else {
				hi = mid;
			}
		}
		return lo;
	}

	/**
	 * Returns the number of indexed entries, which exceeds the number of assemblies: an assembly is
	 * indexed under both its taxon and its species and under both size columns.
	 *
	 * @return the number of indexed entries
	 */
	public int size() {
		return taxids.length;
	}
}
