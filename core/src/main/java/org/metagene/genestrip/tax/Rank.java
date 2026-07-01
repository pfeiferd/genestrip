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
package org.metagene.genestrip.tax;

import org.metagene.genestrip.util.ByteArrayUtil;

/**
 * Taxonomic ranks used in the taxonomy tree, ordered from higher to lower ranks.
 * Besides the standard NCBI ranks it defines a few artificial ranks
 * (REFINED, DATA, FILE, ID) used to track where k-mers originate from. Each rank
 * carries a numeric {@code level} that encodes its relative position independently
 * of the enum ordinal, so that new ranks can be added without breaking already
 * stored databases.
 */
public enum Rank {
	/** The cellular root rank. */
	CELLULAR_ROOT("cellular root"),
	/** The acellular root rank. */
	ACELLULAR_ROOT("acellular root"),
	/** The superkingdom rank. */
	SUPERKINGDOM("superkingdom"),
	/** The domain rank. */
	DOMAIN("domain"),
	/** The realm rank. */
	REALM("realm"),
	/** The kingdom rank. */
	KINGDOM("kingdom"),
	/** The phylum rank. */
	PHYLUM("phylum"),
	/** The subphylum rank. */
	SUBPHYLUM("subphylum"),
	/** The superclass rank. */
	SUPERCLASS("superclass"),
	/** The class rank. */
	CLASS("class"),
	/** The subclass rank. */
	SUBCLASS("subclass"),
	/** The superorder rank. */
	SUPERORDER("superorder"),
	/** The order rank. */
	ORDER("order"),
	/** The suborder rank. */
	SUBORDER("suborder"),
	/** The superfamily rank. */
	SUPERFAMILY("superfamily"),
	/** The family rank. */
	FAMILY("family"),
	/** The subfamily rank. */
	SUBFAMILY("subfamily"),
	/** The tribe rank. */
	TRIBE("tribe"),
	/** The genus rank. */
	GENUS("genus"),
	/** The subgenus rank. */
	SUBGENUS("subgenus"),
	/** The species group rank. */
	SPECIES_GROUP("species group"),
	/** The species rank. */
	SPECIES("species"),
	/** The varietas rank. */
	VARIETAS("varietas"),
	/** The subspecies rank. */
	SUBSPECIES("subspecies"),
	/** The serogroup rank. */
	SEROGROUP("serogroup"),
	/** The biotype rank. */
	BIOTYPE("biotype"),
	/** The strain rank. */
	STRAIN("strain"),
	/** The serotype rank. */
	SEROTYPE("serotype"),
	/** The genotype rank. */
	GENOTYPE("genotype"),
	/** The forma rank. */
	FORMA("forma"),
	/** The forma specialis rank. */
	FORMA_SPECIALIS("forma specialis"),
	/** The isolate rank. */
	ISOLATE("isolate"),
	// "clade" is not a specific rank it is just a collection of related taxons (could be on any level apparently).
	/** The clade rank; a collection of related taxons with an indeterminate level. */
	CLADE("clade", -1),
	// "no rank" is sort of a wild card in the taxonomic tree for both intermediate nodes or leaf notes - typically under in lower ranks.
	/** The wildcard "no rank" rank, with an indeterminate level. */
	NO_RANK("no rank", -1),
	// Recently (newly found ranks) to be added here to preserve order via ordinal from further up.
	/** The subkingdom rank, ordered just below {@link #KINGDOM}. */
	SUBKINGDOM("subkingdom", KINGDOM.level + 10),
	/** The section rank, ordered just below {@link #GENUS}. */
	SECTION("section", GENUS.level + 10),
	// Additional, artificial ranks for tracking where k-mers originate from:
	// They are always the second lowest or the lowest in the taxonomy.
	/** Artificial rank marking refined k-mer origin, near the bottom of the taxonomy. */
	REFINED("REFINED"),
	/** Artificial rank marking a data source of k-mer origin, near the bottom of the taxonomy. */
	DATA("DATA"),
	/** Artificial rank marking the file of k-mer origin, near the bottom of the taxonomy. */
	FILE("FILE"),
	/** Artificial rank marking the id of k-mer origin, at the bottom of the taxonomy. */
	ID("ID");

	/** Level value indicating that a rank has no well-defined level in the taxonomy. */
	public static final int INDETERMINATE_LEVEL = -1;

	private static final Rank[] VALUES = Rank.values();

	private final String name;
	// The level allows for adding new ranks without messing up then rank encoding via ordinals
	// for already stored databases.
	private final int level;

	/**
	 * Whether this rank has no well-defined level (such as {@code CLADE} or
	 * {@code NO_RANK}) and hence cannot be ordered against other ranks.
	 *
	 * @return {@code true} if this rank has no well-defined level
	 */
	public boolean isIndeterminate() {
		return level == INDETERMINATE_LEVEL;
	}

	/**
	 * Returns the rank whose name equals the given string, or {@code null} if none matches.
	 *
	 * @param name the rank name to look up
	 * @return the matching rank, or {@code null} if none matches
	 */
	public static Rank byName(String name) {
		for (Rank r : VALUES) {
			if (r.name.equals(name)) {
				return r;
			}
		}
		return null;
	}
	
	/**
	 * Returns the rank whose name matches the given byte sub-array range, or {@code null} if none matches.
	 *
	 * @param outerArray the byte array containing the rank name
	 * @param start the start index (inclusive) of the name within the array
	 * @param end the end index (exclusive) of the name within the array
	 * @return the matching rank, or {@code null} if none matches
	 */
	public static Rank getRankFromBytes(byte[] outerArray, int start, int end) {
		for (Rank rank : Rank.VALUES) {
			if (ByteArrayUtil.equals(outerArray, start, end, rank.getName())) {
				return rank;
			}
		}
		return null;
	}
	
	private Rank(String name) {
		this.name = name;
		this.level = ordinal() * 20;
	}

	private Rank(String name, int level) {
		this.name = name;
		this.level = level;
	}

	/**
	 * Returns this rank's numeric level.
	 *
	 * @return this rank's numeric level
	 */
	public int getLevel() {
		return level;
	}

	/**
	 * Returns this rank's name.
	 *
	 * @return this rank's name
	 */
	public String getName() {
		return name;
	}
	
	@Override
	public String toString() {
		return name;
	}


	/**
	 * Returns the rank for the given enum ordinal, or {@code null} if {@code i} is -1.
	 *
	 * @param i the enum ordinal, or -1 for none
	 * @return the rank for the given ordinal, or {@code null} if {@code i} is -1
	 */
	public static Rank byOrdinal(int i) {
		return i == -1 ? null : VALUES[i];
	}

	/**
	 * Whether this rank is strictly below (lower/more specific than) the given rank.
	 * Returns {@code false} if the ranks are not comparable.
	 *
	 * @param rank the rank to compare against
	 * @return {@code true} if this rank is strictly below the given rank
	 */
	public boolean isBelow(Rank rank) {
		if (!isComparableTo(rank)) {
			return false;
		}
		return this.level > rank.level;
	}

	/**
	 * Whether this rank is strictly above (higher/less specific than) the given rank.
	 * Returns {@code false} if the ranks are not comparable.
	 *
	 * @param rank the rank to compare against
	 * @return {@code true} if this rank is strictly above the given rank
	 */
	public boolean isAbove(Rank rank) {
		if (!isComparableTo(rank)) {
			return false;
		}
		return this.level < rank.level;
	}

	/**
	 * Whether this rank and the given rank both have a defined level and can therefore be ordered.
	 *
	 * @param rank the rank to compare against
	 * @return {@code true} if both ranks have a defined level and can be ordered
	 */
	public boolean isComparableTo(Rank rank) {
		return !isIndeterminate() && rank != null && !rank.isIndeterminate();
	}
}