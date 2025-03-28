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
	SUPERKINGDOM("superkingdom"), KINGDOM("kingdom"), PHYLUM("phylum"), SUBPHYLUM("subphylum"),
	SUPERCLASS("superclass"), CLASS("class"), SUBCLASS("subclass"), SUPERORDER("superorder"), ORDER("order"),
	SUBORDER("suborder"), SUPERFAMILY("superfamily"), FAMILY("family"), SUBFAMILY("subfamily"),
	TRIBE("tribe"), GENUS("genus"), SUBGENUS("subgenus"), SPECIES_GROUP("species group"), SPECIES("species"), VARIETAS("varietas"),
	SUBSPECIES("subspecies"), SEROGROUP("serogroup"), BIOTYPE("biotype"), STRAIN("strain"), SEROTYPE("serotype"),
	GENOTYPE("genotype"), FORMA("forma"), FORMA_SPECIALIS("forma specialis"), ISOLATE("isolate"),
	NO_RANK("no rank");

	private static Rank[] VALUES = Rank.values();

	private String name;

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
		if (rank == null || NO_RANK.equals(this) || NO_RANK.equals(rank)) {
			return false;
		}
		return this.ordinal() > rank.ordinal();
	}	
}