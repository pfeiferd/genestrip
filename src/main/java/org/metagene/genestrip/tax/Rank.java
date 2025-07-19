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

public enum Rank {
	CELLULAR_ROOT("cellular root"), ACELLULAR_ROOT("acellular root"), SUPERKINGDOM("superkingdom"), DOMAIN("domain"), REALM("realm"), KINGDOM("kingdom"), PHYLUM("phylum"), SUBPHYLUM("subphylum"),
	SUPERCLASS("superclass"), CLASS("class"), SUBCLASS("subclass"), SUPERORDER("superorder"), ORDER("order"),
	SUBORDER("suborder"), SUPERFAMILY("superfamily"), FAMILY("family"), SUBFAMILY("subfamily"),
	TRIBE("tribe"), GENUS("genus"), SUBGENUS("subgenus"), SPECIES_GROUP("species group"), SPECIES("species"), VARIETAS("varietas"),
	SUBSPECIES("subspecies"), SEROGROUP("serogroup"), BIOTYPE("biotype"), STRAIN("strain"), SEROTYPE("serotype"),
	GENOTYPE("genotype"), FORMA("forma"), FORMA_SPECIALIS("forma specialis"), ISOLATE("isolate"),
	// "clade" is not a specific rank it is just a collection of related taxons (could be on any level apparently).
	CLADE("clade", -1),
	// "no rank" is sort of a wild card in the taxonomic tree for both intermediate nodes or leaf notes - typically under in lower ranks.
	NO_RANK("no rank", -1),
	// Recently (newly found ranks) to be added here to preserve order via ordinal from further up.
	SUBKINGDOM("subkingdom", KINGDOM.ordinal() + 10),
	SECTION("section", GENUS.ordinal() + 10);

	public static final int INDETERMINATE_LEVEL = -1;

	private static Rank[] VALUES = Rank.values();

	private final String name;
	// The level allows for adding new ranks without messing up then rank encoding via ordinals
	// for already stored databases.
	private final int level;

	public boolean isIndeterminate() {
		return level == INDETERMINATE_LEVEL;
	}

	public static Rank byName(String name) {
		for (Rank r : VALUES) {
			if (r.name.equals(name)) {
				return r;
			}
		}
		return null;
	}
	
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

	public int getLevel() {
		return level;
	}

	public String getName() {
		return name;
	}
	
	@Override
	public String toString() {
		return name;
	}


	public static Rank byOrdinal(int i) {
		return i == -1 ? null : VALUES[i];
	}

	public boolean isBelow(Rank rank) {
		if (!isComparableTo(rank)) {
			return false;
		}
		return this.level > rank.level;
	}

	public boolean isAbove(Rank rank) {
		if (!isComparableTo(rank)) {
			return false;
		}
		return this.level < rank.level;
	}

	public boolean isComparableTo(Rank rank) {
		return !isIndeterminate() && rank != null && !rank.isIndeterminate();
	}
}