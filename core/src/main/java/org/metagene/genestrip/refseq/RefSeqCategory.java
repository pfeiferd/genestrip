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

import org.metagene.genestrip.util.ByteArrayUtil;

/**
 * The NCBI RefSeq release categories (each corresponding to a directory) that genomes are organized
 * into.
 */
public enum RefSeqCategory {
	/** The archaea RefSeq category. */
	ARCHAEA("archaea", true),
	/** The bacteria RefSeq category. */
	BACTERIA("bacteria", true),
	/** The whole release in one directory: every other category's sequences, not complete genomes. */
	COMPLETE("complete", false),
	/** The fungi RefSeq category. */
	FUNGI("fungi", false),
	/** The invertebrate RefSeq category. */
	INVERTEBRATE("invertebrate", false),
	/** The mitochondrion RefSeq category. */
	MITOCHONDRION("mitochondrion", true),
	/** The other RefSeq category. */
	OTHER("other", true),
	/** The plant RefSeq category. */
	PLANT("plant", false),
	/** The plasmid RefSeq category. */
	PLASMID("plasmid", true),
	/** The plastid RefSeq category. */
	PLASTID("plastid", true),
	/** The protozoa RefSeq category. */
	PROTOZOA("protozoa", false),
	/** The mammalian vertebrate RefSeq category. */
	VERTEBRATE_MAMMALIAN("vertebrate_mammalian", false),
	/** The other vertebrate RefSeq category. */
	VERTEBRATE_OTHER("vertebrate_other", false),
	/** The viral RefSeq category. */
	VIRAL("viral", true);

	/** The RefSeq directory name associated with this category. */
	private String directory;

	/** Whether the genome of a sequence in this category can be told from its accession. */
	private boolean genomeKeyed;

	/**
	 * Creates a category for the given RefSeq directory name.
	 *
	 * @param directory   the RefSeq directory name.
	 * @param genomeKeyed whether the genome of a sequence can be told from its accession.
	 */
	RefSeqCategory(String directory, boolean genomeKeyed) {
		this.directory = directory;
		this.genomeKeyed = genomeKeyed;
	}

	/**
	 * Returns whether {@link GenomeKeyTrie#genomeKeyLength(byte[], int, int)} finds the genome of a
	 * sequence in this category from its accession, which is what {@code genomeNodes} needs.
	 * <p>
	 * It does where a genome is a WGS project or one replicon: prokaryotes, organelles, plasmids and
	 * viruses. It does not for eukaryotes. Their RefSeq genomes are chromosomes ({@code NC_}) and
	 * scaffolds ({@code NW_}) with consecutive accessions and no shared prefix, so every one of them
	 * becomes a genome of its own: in release 236 the protozoa category has 191,236 genomic sequences,
	 * and as many keys, for 129 assemblies. The catalog names no assembly, so nothing in the accession
	 * can mend this. {@link #COMPLETE} holds the eukaryotes too.
	 *
	 * @return whether the accession identifies the genome in this category.
	 */
	public boolean isGenomeKeyed() {
		return genomeKeyed;
	}

	/**
	 * Returns whether {@code genomeNodes=true} takes all sequences of one taxon in this category's
	 * release files for one genome, instead of telling the genomes apart by accession. That is so for
	 * every category that is not {@link #isGenomeKeyed()}, except {@link #COMPLETE}: it holds the
	 * prokaryotes as well, so it is decided there per taxon, by the catalog - see
	 * {@link #isOneGenomePerTaxonRow(byte[], int, int)}.
	 * <p>
	 * It holds because RefSeq carries one assembly per taxon for the eukaryotes: in release 236 the
	 * assembly summary lists two for 7 of their 2,419 taxa and one for all others.
	 *
	 * @return whether this category's release files get one genome per taxon.
	 */
	public boolean isOneGenomePerTaxon() {
		return !genomeKeyed && this != COMPLETE;
	}

	/**
	 * Returns whether a row of the accession catalog belongs to a category that is
	 * {@link #isOneGenomePerTaxon()}, from the row's release directories, such as
	 * {@code complete|protozoa}. Every row lists {@code complete} beside the category proper, and every
	 * row of a eukaryote lists its eukaryotic category, organelles included
	 * ({@code complete|mitochondrion|protozoa}), so this tells which taxa of the {@link #COMPLETE}
	 * release files have one genome per taxon.
	 *
	 * @param row   the buffer holding the catalog row
	 * @param start the start index of the directories column (inclusive)
	 * @param end   the end index of the directories column (exclusive)
	 * @return whether the row belongs to a category with one genome per taxon
	 */
	public static boolean isOneGenomePerTaxonRow(byte[] row, int start, int end) {
		for (RefSeqCategory cat : values()) {
			if (cat.isOneGenomePerTaxon() && ByteArrayUtil.indexOf(row, start, end, cat.getDirectory()) != -1) {
				return true;
			}
		}
		return false;
	}

	/**
	 * Returns the category of a RefSeq release file from its name, such as
	 * {@code protozoa.1.1.genomic.fna.gz}, or null if the name starts with no category.
	 *
	 * @param fileName the name of the release file.
	 * @return the matching category, or null if none matches.
	 */
	public static RefSeqCategory fromFileName(String fileName) {
		int pos = fileName.indexOf('.');
		return pos < 0 ? null : fromDirectoryString(fileName.substring(0, pos));
	}

	/**
	 * Returns the RefSeq directory name associated with this category.
	 *
	 * @return the directory name.
	 */
	public String getDirectory() {
		return directory;
	}
	
	/**
	 * Returns the category whose directory name equals the given string, or null if none matches.
	 *
	 * @param category the RefSeq directory name to look up.
	 * @return the matching category, or null if none matches.
	 */
	public static RefSeqCategory fromDirectoryString(String category) {
		for (RefSeqCategory cat : RefSeqCategory.values()) {
			if (cat.getDirectory().equals(category)) {
				return cat;
			}
		}
		return null;
	}	
}
