/*
 * Genestrip
 */
package org.metagene.genestrip.refseq;

import org.metagene.genestrip.genbank.AssemblySummaryReader.AssemblyQuality;

/**
 * What is known about the assembly a sequence belongs to.
 * <p>
 * Two things are kept and they answer different questions. The key is the one every sequence of an
 * assembly is filed under, so that a chromosome and its plasmids count as one genome rather than as
 * three. The level says how finished that assembly is, which decides whether the genome may carry a
 * topology: a tree inferred from genomes that are short by an unknown amount has the shortness in it,
 * whereas one inferred from complete genomes and then extended by placing the rest does not.
 * <p>
 * The level is carried here rather than consumed and discarded because filtering and annotating are
 * different uses. A build that restricts itself to finished genomes needs only to know whether an
 * assembly was admitted; a build that keeps everything and wants to know which of its genomes are
 * complete needs the answer itself, and the same pass produces both.
 */
public class AssemblyInfo {

	private final byte[] key;

	private final AssemblyQuality level;

	/**
	 * Creates the record.
	 *
	 * @param key the genome key every sequence of the assembly is filed under
	 * @param level the assembly level the summary records
	 */
	public AssemblyInfo(byte[] key, AssemblyQuality level) {
		this.key = key;
		this.level = level;
	}

	/**
	 * Returns the genome key shared by every sequence of this assembly.
	 *
	 * @return the assembly's genome key
	 */
	public byte[] getKey() {
		return key;
	}

	/**
	 * Returns the assembly level as the summary records it.
	 *
	 * @return the assembly level
	 */
	public AssemblyQuality getLevel() {
		return level;
	}

	/**
	 * Returns whether the assembly is complete, i.e. every replicon closed end to end. This is the
	 * question a backbone asks: a complete assembly holds the whole of its organism's genome, so its
	 * distances are not inflated by anything missing and it may be used to infer a topology.
	 *
	 * @return whether the assembly is complete
	 */
	public boolean isComplete() {
		return AssemblySizeIndex.COMPLETE_ONLY.contains(level);
	}

	/**
	 * Returns whether the assembly is chromosome-level rather than complete.
	 * <p>
	 * Such an assembly falls about one per cent short of the complete genomes of its own species. That
	 * per cent is not recorded anywhere: the summary's two size columns count only the gaps the
	 * assembler declared, whose median is a fiftieth of the shortfall a same-species comparison
	 * implies. So this is a category and not a measurement, and the shortfall it stands for should be
	 * displayed rather than corrected for.
	 *
	 * @return whether the assembly is chromosome-level
	 */
	public boolean isChromosome() {
		return AssemblySizeIndex.CHROMOSOME_LEVEL.contains(level);
	}
}
