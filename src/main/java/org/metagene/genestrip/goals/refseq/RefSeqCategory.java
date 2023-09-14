package org.metagene.genestrip.goals.refseq;

public enum RefSeqCategory {
	archaea("archaea"),
	bacteria("bacteria"),
	complete("complete"),
	fungi("fungi"),
	invertebrate("invertebrate"),
	mitochondrion("mitochondrion"),
	other("other"),
	plant("plant"),
	plasmid("plasmid"),
	plastid("plastid"),
	protozoa("protozoa"),
	vertebrate_mammalian("vertebrate_mammalian"),
	vertebrate_other("vertebrate_other"),
	viral("viral"); 

	private String directory;
	
	RefSeqCategory(String directory) {
		this.directory = directory;
	}
	
	public String getDirectory() {
		return directory;
	}
}
