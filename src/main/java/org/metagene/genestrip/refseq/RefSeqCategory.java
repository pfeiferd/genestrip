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

public enum RefSeqCategory {
	ARCHAEA("archaea"), BACTERIA("bacteria"), COMPLETE("complete"), FUNGI("fungi"), INVERTEBRATE("invertebrate"),
	MITOCHONDRION("mitochondrion"), OTHER("other"), PLANT("plant"), PLASMID("plasmid"), PLASTID("plastid"),
	PROTOZOA("protozoa"), VERTEBRATE_MAMMALIAN("vertebrate_mammalian"), VERTEBRATE_OTHER("vertebrate_other"),
	VIRAL("viral");

	private String directory;

	RefSeqCategory(String directory) {
		this.directory = directory;
	}

	public String getDirectory() {
		return directory;
	}
	
	public static RefSeqCategory fromDiregoryString(String category) {
		for (RefSeqCategory cat : RefSeqCategory.values()) {
			if (cat.getDirectory().equals(category)) {
				return cat;
			}
		}
		return null;
	}	
}