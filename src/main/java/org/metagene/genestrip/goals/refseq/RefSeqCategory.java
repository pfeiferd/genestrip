package org.metagene.genestrip.goals.refseq;

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
