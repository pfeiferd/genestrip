package org.metagene.genestrip;

import java.io.File;

import org.metagene.genestrip.make.FileDownloadGoal.DownloadProject;
import org.metagene.genestrip.tax.AssemblySummaryReader.FTPEntryQuality;

public class GSProject implements DownloadProject {
	private final GSConfig config;
	private final String name;
	private final String krakenDB;
	private final FTPEntryQuality fastaQuality;
	private final int kMserSize;
	private final File fastqFile;
	private final File resFolder;

	public GSProject(GSConfig config, String name, FTPEntryQuality fastaQuality, int kMerSize, String krakenDB,
			File fastqFile, File resFolder) {
		this.config = config;
		this.name = name;
		this.fastaQuality = fastaQuality;
		this.kMserSize = kMerSize;
		this.krakenDB = krakenDB;
		this.fastqFile = fastqFile;
		this.resFolder = resFolder != null ? resFolder : new File(getProjectsDir(), name + "/res");
	}

	public File getFastqFile() {
		return fastqFile;
	}
	
	@Override
	public boolean isIgnoreMissingFiles() {
		return getConfig().isIgnoreMissingFastas();
	}
	
	@Override
	public String getBaseFTPURL() {
		return getConfig().getFtpBaseURL();
	}
	
	@Override
	public String getHttpBaseURL() {
		return getConfig().getHttpBaseURL();
	}
	
	@Override
	public boolean isUseHttp() {
		return getConfig().isUseHttp();
	}

	public File getFilteredFastqFile() {
		return new File(getResultsDir(), "Filtered_" + fastqFile.getName());
	}

	public File getDumpFastqFile() {
		return getConfig().isWriteDumpedFastq() ? new File(getResultsDir(), "Dumped_" + fastqFile.getName()) : null;
	}

	public File getTaxCountsFile() {
		return new File(getResultsDir(), "Counts_" + fastqFile.getName() + ".csv");
	}

	public GSConfig getConfig() {
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
		return resFolder;
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
