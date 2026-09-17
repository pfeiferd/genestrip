/*
 * Genestrip
 */
package org.metagene.genestrip.goals.refseq;

import java.io.File;
import java.io.IOException;

import org.metagene.genestrip.GSConfigKey;
import org.metagene.genestrip.GSGoalKey;
import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.genbank.AssemblySummaryReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.List;
import java.util.Set;

import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.refseq.AccessionFileProcessor;
import org.metagene.genestrip.refseq.AssemblySizeIndex;
import org.metagene.genestrip.refseq.GenomeKeyTrie;
import org.metagene.genestrip.refseq.AccessionTrie;
import org.metagene.genestrip.genbank.AssemblySummaryReader.AssemblyQuality;
import org.metagene.genestrip.refseq.AssemblyInfo;
import org.metagene.genestrip.GSConfigKey.SeqType;
import org.metagene.genestrip.io.StreamingFileResource;
import org.metagene.genestrip.io.StreamingResource;
import org.metagene.genestrip.refseq.RefSeqCategory;

/**
 * Provides what NCBI's assembly summary knows about the assemblies a database is built from.
 * <p>
 * The accession catalog names sequences and carries no assembly accession, so nothing in it says
 * whether a sequence belongs to a finished genome or to a draft. The summary says so, per assembly,
 * in its {@code assembly_level} column. {@link AssemblySizeIndex} joins the two on the taxon and the
 * number of bases, which is the one thing they share, and this goal makes that index available.
 * <p>
 * It is an object goal and so is made only when something asks for it: a build that leaves
 * {@code refseq.genomesOnly} at {@code off} never reads the summary and pays nothing.
 * <p>
 * What the goal offers is the grouping: for every sequence of a finished assembly, the key of that
 * assembly, so that a caller can both tell which sequences to keep and count assemblies rather than
 * the replicons they consist of.
 * <p>
 * Turning an accession into an assembly is what the catalog pass does: it is ordered by taxon and,
 * within a taxon, by accession, so a run of consecutive accessions can be gathered and summed in a
 * single streaming pass with a window of a handful of records and no map, and the sum looked up in
 * the index.
 * <p>
 * What the goal does not offer is the rest of the summary's columns, because nothing consumes them.
 * Keeping them belongs here, beside the grouping, and the index already reads every level when it is
 * built with {@link AssemblySizeIndex#ALL_LEVELS}, which is the other half of what such a consumer would need.
 *
 * @param <P> the type of the project this goal belongs to.
 */
public class AssemblyMetadataGoal<P extends GSProject> extends ObjectGoal<AccessionTrie<AssemblyInfo>, P> {

	/** The longest run of consecutive accessions taken to be one assembly. */
	private static final int MAX_RUN = 8;

			/**
	 * One catalog entry, held while its run is still open.
	 * <p>
	 * The buffer is reused: a run is at most {@link #MAX_RUN} entries, so the slots are made
	 * once and refilled. The catalog has hundreds of millions of rows and every one of them
	 * passes through here, which is the same reason it is read into one line buffer rather
	 * than a string per line.
	 */
	private static class Entry {
		final byte[] accession = new byte[MAX_ACCESSION];
		int length0;
		long number;
		long length;

		void set(byte[] target, int start, int end, long number, long length) {
			this.length0 = Math.min(end - start, MAX_ACCESSION);
			System.arraycopy(target, start, accession, 0, length0);
			this.number = number;
			this.length = length;
		}
	}

	/** Longest accession held; RefSeq accessions are far shorter than this. */
	private static final int MAX_ACCESSION = 64;

	private final RefSeqCatalogDownloadGoal catalogGoal;
	private final ObjectGoal<Set<RefSeqCategory>, P> categoriesGoal;

	/**
	 * Creates the goal, wiring the catalog and categories goals it reads.
	 *
	 * @param project the project
	 * @param categoriesGoal goal providing the selected RefSeq categories
	 * @param catalogGoal goal providing the downloaded RefSeq catalog
	 * @param deps additional goal dependencies
	 */
	@SafeVarargs
	public AssemblyMetadataGoal(P project, ObjectGoal<Set<RefSeqCategory>, P> categoriesGoal,
			RefSeqCatalogDownloadGoal catalogGoal, Goal<P>... deps) {
		super(project, GSGoalKey.ASSEMBLYMETA, Goal.append(deps, categoriesGoal, catalogGoal));
		this.categoriesGoal = categoriesGoal;
		this.catalogGoal = catalogGoal;
	}

	/**
	 * Returns how many bases a sum of sequence lengths may differ from an assembly's recorded length
	 * and still be taken for that assembly.
	 *
	 * @return the configured gap
	 */
	public int getGap() {
		return intConfigValue(GSConfigKey.ASSEMBLY_METADATA_GAP);
	}

	/**
	 * Returns what is known about the assembly the given accession belongs to, or {@code null} where
	 * none was matched. The goal is made if it has not been already, so this may be called without
	 * calling {@link #get()} first.
	 * <p>
	 * The accession is given as it appears in the catalog or in a FASTA header, version suffix and
	 * all -- {@code NC_001911.1}. Deriving the key the assemblies are filed under is this method's
	 * job, and it is not a trivial one: the version is cut, and a whole-genome shotgun accession is
	 * cut back to its project prefix so that the contigs of one draft share a key.
	 * <p>
	 * A {@code null} answer says that no admitted assembly was matched, which is <em>not</em> the same
	 * as saying the genome is incomplete. It happens for a draft, for a genome whose sum of sequence
	 * lengths did not meet an assembly within {@code refseq.assemblyMetadataGap}, and for an accession
	 * that is not genomic at all. A caller selecting the genomes that may carry a topology should
	 * treat {@code null} as "not established as complete" and not as "established as not complete".
	 *
	 * @param target the buffer holding the accession
	 * @param start the start offset of the accession
	 * @param end the end offset of the accession
	 * @return what is known about the assembly, or {@code null}
	 */
	public AssemblyInfo getAssemblyInfo(byte[] target, int start, int end) {
		int keyEnd = start + GenomeKeyTrie.genomeKeyLength(target, start, end);
		return get().get(target, start, keyEnd);
	}

	/**
	 * Returns what is known about the assembly the given accession belongs to, or {@code null} where
	 * none was matched. See {@link #getAssemblyInfo(byte[], int, int)} for what {@code null} means.
	 * <p>
	 * This makes a byte array per call and is meant for asking about one accession. Code walking a
	 * whole database should use the byte-array form against a buffer it already has, as the goals in
	 * this package do.
	 *
	 * @param accession the accession, e.g. {@code NC_001911.1}
	 * @return what is known about the assembly, or {@code null}
	 */
	public AssemblyInfo getAssemblyInfo(String accession) {
		byte[] target = new byte[accession.length()];
		for (int i = 0; i < target.length; i++) {
			target[i] = (byte) accession.charAt(i);
		}
		return getAssemblyInfo(target, 0, target.length);
	}

	@Override
	protected void doMakeThis() {
		File dir = getProject().getCommon().getRefSeqDir();
		File summary = new File(dir, AssemblySummaryReader.ASSEMLY_SUM_REFSEQ);
		if (!summary.exists()) {
			throw new IllegalStateException(
					"The assembly summary " + summary + " is needed for assembly metadata but is not there.");
		}
		final AssemblySizeIndex index;
		try {
			// The index is built over exactly the levels the setting admits, so that everything in
			// it is something to keep and a search needs no screening afterwards. Screening after
			// the search would lose a run whose admitted assembly lay within the gap but which an
			// unadmitted one happened to lie nearer to.
			// Complete assemblies are always indexed, because whether a genome is complete is what a
			// consumer of this goal asks even when nothing is being filtered. Chromosome-level ones
			// are added only where the filter admits them, since indexing a level the build will not
			// keep would let it win a search from the level the build does keep.
			GSConfigKey.GenomesOnly genomesOnly = (GSConfigKey.GenomesOnly) configValue(GSConfigKey.GENOMES_ONLY);
			index = new AssemblySizeIndex(summary,
					genomesOnly == GSConfigKey.GenomesOnly.CHROMOSOME
							? AssemblySizeIndex.COMPLETE_OR_CHROMOSOME
							: AssemblySizeIndex.COMPLETE_ONLY,
					genomesOnly == GSConfigKey.GenomesOnly.PREF_REF);
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
		// An accession-keyed trie, not a plain DigitTrie: the keys are accessions, and a DigitTrie
		// maps every letter outside its range, so a write would find no node and fail.
		final AccessionTrie<AssemblyInfo> complete = new AccessionTrie<AssemblyInfo>();
		final int gap = getGap();

		AccessionFileProcessor processor = new AccessionFileProcessor(categoriesGoal.get(),
				(SeqType) configValue(GSConfigKey.SEQ_TYPE),
				(List<GSConfigKey.RefSeqStatus>) configValue(GSConfigKey.RES_SEQ_STATUS)) {

			private final Entry[] run = newRun();
			private int runSize;
			/**
			 * The assemblies already accounted for. An assembly is one genome however many runs
			 * happen to match it, and several do: a genome's chromosome alone can match another
			 * assembly of the species more closely than the whole genome matches its own, and in a
			 * densely sequenced species the recorded lengths lie close enough together that a sum
			 * meets one by chance. Counting each sighting would count one genome many times.
			 */
			private final BitSet seenAssemblies = new BitSet(index.getAssemblyCount());
			/** How many assemblies were recognised, and how many sequences they account for. */
			private int assemblies;
			private int sequences;
			/** Genomic rows seen, so that the share recognised can be told. */
			private int genomicRows;
			private int noLength, noTaxId, noDigits, runsFlushed, lookups, hits;
			private long maxRun;
			private int runTaxId = -1;
			/** The letters of the run's accessions, held as bytes so that no string is made per row. */
			private final byte[] runPrefix = new byte[MAX_ACCESSION];
			private int runPrefixLength;

			private Entry[] newRun() {
				Entry[] r = new Entry[MAX_RUN];
				for (int i = 0; i < r.length; i++) {
					r[i] = new Entry();
				}
				return r;
			}

			@Override
			protected void handleEntry(byte[] target, int taxIdEnd, int accessionStart, int accessionEnd) {
				// Never called: the overload below is, and it does not delegate.
			}

			@Override
			protected void handleEntry(byte[] target, int taxIdEnd, int accessionStart, int accessionEnd,
					int lengthStart, int lengthEnd) {
				genomicRows++;
				long length = parseLength(target, lengthStart, lengthEnd);
				if (length <= 0) {
					noLength++;
					return;
				}
				int taxId = (int) parseLength(target, 0, taxIdEnd);
				if (taxId <= 0) {
					noTaxId++;
					return;
				}
				// The accession splits into letters and digits after the two-letter type and underscore.
				// Sequences of one assembly are numbered consecutively within one such letter group, and
				// the catalog lists them in that order, so a run is recognised as it streams past.
				int p = accessionStart + 3;
				int letters = p;
				while (letters < accessionEnd && target[letters] >= 'A' && target[letters] <= 'Z') {
					letters++;
				}
				int digits = letters;
				while (digits < accessionEnd && target[digits] >= '0' && target[digits] <= '9') {
					digits++;
				}
				if (digits == letters) {
					noDigits++;
					return;
				}
				long number = parseLength(target, letters, digits);
				if (number < 0) {
					return;
				}
				int prefixLength = letters - accessionStart;
				boolean continues = taxId == runTaxId && runSize > 0
						&& samePrefix(target, accessionStart, prefixLength)
						&& number == run[runSize - 1].number + 1;
				if (!continues) {
					flush();
					runTaxId = taxId;
					setPrefix(target, accessionStart, prefixLength);
				}
				run[runSize++].set(target, accessionStart, accessionEnd, number, length);
				if (runSize == MAX_RUN) {
					flush();
					runTaxId = taxId;
					setPrefix(target, accessionStart, prefixLength);
				}
			}

			/**
			 * Decides the open run and admits whatever of it belongs to a finished genome.
			 * <p>
			 * The sequences of one assembly sum to its recorded length, so the run is walked from the
			 * front and, at each position, the grouping whose sum comes nearest a recorded length wins.
			 * A grouping that turns out to be a finished assembly has the genome keys of its members
			 * admitted, which for a replicon is the accession itself and for a shotgun project is the
			 * project, so admitting a project's master admits its contigs with it.
			 */
			/** Whether the accession's letters are the ones the open run carries. */
			private boolean samePrefix(byte[] target, int start, int length) {
				if (length != runPrefixLength) {
					return false;
				}
				for (int i = 0; i < length; i++) {
					if (target[start + i] != runPrefix[i]) {
						return false;
					}
				}
				return true;
			}

			private void setPrefix(byte[] target, int start, int length) {
				runPrefixLength = Math.min(length, MAX_ACCESSION);
				System.arraycopy(target, start, runPrefix, 0, runPrefixLength);
			}

			private void flush() {
				int n = runSize;
				if (n > 0) {
					runsFlushed++;
					maxRun = Math.max(maxRun, n);
				}
				int i = 0;
				while (i < n) {
					int bestLen = 0;
					long bestDist = Long.MAX_VALUE;
					int bestId = -1;
					AssemblyQuality bestLevel = null;
					boolean bestRef = false;
					long sum = 0;
					for (int len = 1; len <= n - i; len++) {
						Entry e = run[i + len - 1];
						if (e.number != run[i].number + len - 1) {
							break;
						}
						sum += e.length;
						// The index holds admitted assemblies only, so a match is one of those
						// and the nearest is the best grouping of this run. Screening the level
						// afterwards instead would lose a run whose admitted assembly was in range
						// but a draft happened to lie nearer.
						lookups++;
						AssemblySizeIndex.Match m = index.nearest(runTaxId, sum, gap);
						if (m != null) {
							hits++;
						}
						if (m != null && m.getDistance() < bestDist && !seenAssemblies.get(m.getAssemblyId())) {
							bestDist = m.getDistance();
							bestId = m.getAssemblyId();
							bestLevel = m.getLevel();
							bestRef = m.isReference();
							bestLen = len;
						}
					}
					if (bestLen == 0) {
						i++;
						continue;
					}
					seenAssemblies.set(bestId);
					// Every sequence of the assembly is filed under the same key, that of the first
					// of them, so that a chromosome and its plasmids count as one genome rather than
					// as three. Without this the cap of maxGenomesPerTaxid would spend its quota on
					// replicons: a finished genome arrives as several accessions, each of which the
					// genome key leaves standing for itself.
					Entry leader = run[i];
					AssemblyInfo info = new AssemblyInfo(Arrays.copyOf(leader.accession,
							GenomeKeyTrie.genomeKeyLength(leader.accession, 0, leader.length0)), bestLevel, bestRef);
					assemblies++;
					sequences += bestLen;
					for (int j = i; j < i + bestLen; j++) {
						Entry e = run[j];
						int keyEnd = GenomeKeyTrie.genomeKeyLength(e.accession, 0, e.length0);
						complete.set(e.accession, 0, keyEnd, info);
					}
					i += bestLen;
				}
				runSize = 0;
			}

			@Override
			public void processCatalog(StreamingResource catalogFile) {
				super.processCatalog(catalogFile);
				flush();
				if (getLogger().isInfoEnabled()) {
					getLogger().info("Assemblies recognised: " + assemblies);
					getLogger().info("Sequences they account for: " + sequences);
					getLogger().info("Genomic catalog rows seen: " + genomicRows);
					getLogger().info("Dropped - no length: " + noLength + ", no taxon: " + noTaxId
							+ ", no digits: " + noDigits);
					getLogger().info("Runs flushed: " + runsFlushed + ", longest: " + maxRun);
					getLogger().info("Index lookups: " + lookups + ", of which hit: " + hits);
					getLogger().info("Index entries: " + index.size());
				}
			}
		};
		processor.processCatalog(new StreamingFileResource(catalogGoal.getCatalogFile()));
		set(complete);
	}
}
