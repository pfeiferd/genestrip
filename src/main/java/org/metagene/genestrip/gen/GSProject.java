package org.metagene.genestrip.gen;

import java.io.File;

import org.metagene.genestrip.tax.AssemblySummaryReader.FTPEntryQuality;

public class GSProject {
	private final Config config;
	private final String name;
	private final String krakenDB;
	private final FTPEntryQuality fastaQuality;
	private final int kMserSize;
	private final File fastqFile;

	public GSProject(Config config, String name, FTPEntryQuality fastaQuality, int kMerSize, String krakenDB, File fastqFile) {
		this.config = config;
		this.name = name;
		this.fastaQuality = fastaQuality;
		this.kMserSize = kMerSize;
		this.krakenDB = krakenDB;
		this.fastqFile = fastqFile;
	}
	
	public File getFastqFile() {
		return fastqFile;
	}
	
	public File getFilteredFastqFile() {
		return new File(getResultsDir(), "Filtered_" + fastqFile.getName());
	}
	
	public Config getConfig() {
		return config;
	}

	public String getKrakenDB() {
		return krakenDB == null ? config.getKrakenDB() : krakenDB;
	}

	public File getBaseDir() {
		return config.getBaseDir();
	}

	public int getkMserSize() {
		return kMserSize <= 0 ? config.getkMerSize() : kMserSize;
	}

	public FTPEntryQuality getFastaQuality() {
		return fastaQuality == null ? config.getFastaQuality() : fastaQuality;
	}

	public String getName() {
		return name;
	}
	
	public File getProjectsDir() {
		return new File(config.getBaseDir(), "/projects");
	}

	public File getFastasDir() {
		return new File(getProjectsDir(), name + "/fastas");
	}

	public File getFastqsDir() {
		return new File(getProjectsDir(), name + "/fastqs");
	}

	public File getFiltersDir() {
		return new File(getProjectsDir(), name + "/filters");
	}

	public File getTaxIdsFile() {
		return new File(getProjectsDir(), name + "/taxids.txt");
	}
	
	public File getResultsDir() {
		return new File(getProjectsDir(), name + "/res");		
	}

	public File getKmerFastqFile() {
		return new File(getFastqsDir(), getName() + "_k" + getkMserSize() + "_fromFastas.fastq.gz");
	}

	public File getKrakenOutFile() {
		return new File(getFastqsDir(), getName() + "_k" + getkMserSize() + "_fromKraken.out.txt");
	}

	public File getFilteredKmerFastqFile() {
		return new File(getFastqsDir(), getName() + "_k" + getkMserSize() + "_fromKraken.fastq.gz");
	}

	public File getTrieFile() {
		return new File(getFiltersDir(), getName() + "_k" + getkMserSize() + "_Trie.ser.gz");
	}

	public File getBloomFilterFile() {
		return new File(getFiltersDir(), getName() + "_k" + getkMserSize() + "_BloomFilter.ser.gz");
	}
}
