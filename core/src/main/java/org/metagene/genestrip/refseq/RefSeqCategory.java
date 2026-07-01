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

/**
 * The NCBI RefSeq release categories (each corresponding to a directory) that genomes are organized
 * into.
 */
public enum RefSeqCategory {
	/** The archaea RefSeq category. */
	ARCHAEA("archaea"),
	/** The bacteria RefSeq category. */
	BACTERIA("bacteria"),
	/** The complete-genomes RefSeq category. */
	COMPLETE("complete"),
	/** The fungi RefSeq category. */
	FUNGI("fungi"),
	/** The invertebrate RefSeq category. */
	INVERTEBRATE("invertebrate"),
	/** The mitochondrion RefSeq category. */
	MITOCHONDRION("mitochondrion"),
	/** The other RefSeq category. */
	OTHER("other"),
	/** The plant RefSeq category. */
	PLANT("plant"),
	/** The plasmid RefSeq category. */
	PLASMID("plasmid"),
	/** The plastid RefSeq category. */
	PLASTID("plastid"),
	/** The protozoa RefSeq category. */
	PROTOZOA("protozoa"),
	/** The mammalian vertebrate RefSeq category. */
	VERTEBRATE_MAMMALIAN("vertebrate_mammalian"),
	/** The other vertebrate RefSeq category. */
	VERTEBRATE_OTHER("vertebrate_other"),
	/** The viral RefSeq category. */
	VIRAL("viral");

	/** The RefSeq directory name associated with this category. */
	private String directory;

	/**
	 * Creates a category for the given RefSeq directory name.
	 *
	 * @param directory the RefSeq directory name.
	 */
	RefSeqCategory(String directory) {
		this.directory = directory;
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
