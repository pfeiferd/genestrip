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
package org.metagene.genestrip;

import java.io.File;

import org.metagene.genestrip.make.FileDownloadGoal.DownloadProject;
import org.metagene.genestrip.tax.AssemblySummaryReader.FTPEntryQuality;

public class GSProject implements DownloadProject {
	public enum FileType {
		FASTQ(".fastq"), FASTA(".fasta"), CSV(".csv"), KRAKEN_OUT(".out"), SER(".ser");

		private final String suffix;

		private FileType(String suffix) {
			this.suffix = suffix;
		}

		public String getSuffix() {
			return suffix;
		}
	}

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
		this.resFolder = resFolder != null ? resFolder : new File(getProjectsDir(), name + "/csv");
	}

	public File getFastqFile() {
		return fastqFile;
	}

	public File getDirForType(FileType type) {
		switch (type) {
		case FASTQ:
			return getFastqsDir();
		case FASTA:
			return getFastasDir();
		case CSV:
			return getResultsDir();
		case KRAKEN_OUT:
			return getKrakenOutDir();
		case SER:
			return getFiltersDir();
		default:
			throw new IllegalArgumentException("Illegal FileType: " + type);
		}
	}

	public File getOutputFile(String goal, FileType type) {
		return getOutputFile(goal, null, type, true);
	}

//	public File getOutputFileForFastq(String goal, FileType type) {
//		return getOutputFile(goal, getFastqFile(), type);
//	}

	public File getOutputFile(String goal, File baseFile, FileType type) {
		return getOutputFile(goal, baseFile, type, true);
	}

	public File getOutputFile(String goal, File baseFile, FileType type, boolean gzip) {
		String baseName = baseFile == null ? "" : baseFile.getName();
		if (baseName.endsWith(".gz")) {
			baseName = baseName.substring(0, baseName.length() - 3);
		}
		else if (baseName.endsWith(".gzip")) {
			baseName = baseName.substring(0, baseName.length() - 5);
		}
		for (FileType ft : FileType.values()) {
			if (baseName.endsWith(ft.getSuffix())) {
				baseName = baseName.substring(0, ft.getSuffix().length());
			}			
		}		
		if (baseName.startsWith(getName() + "_")) {
			baseName.substring(getName().length() + 1);
		}
		if (!baseName.isEmpty()) {
			baseName = "_" + baseName;
		}
		
		return new File(getDirForType(type), getName() + "_" + goal + baseName + type.getSuffix() + (gzip ? ".gz" : ""));
	}

//	public File getFastqKrakenOutFile() {
//		return new File(getKrakenOutDir(), fastqFile.getName() + getkMserSize() + "_fromKraken.out.txt.gz");
//	}

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

//	public File getFilteredFastqFile() {
//		return new File(getFastqsDir(), "filtered_" + getFileNameWithGzEnding(fastqFile));
//	}

//	public File getDumpFastqFile(File fastq) {
//		return getConfig().isWriteDumpedFastq() ? getOutputFile("dumped", fastq, FileType.FASTQ) : null;
//	}

//	private String getFileNameWithGzEnding(File file) {
//		String name = file.getName();
//		if (name.endsWith(".gz") || name.endsWith(".gzip")) {
//			return name;
//		}
//		return name + ".gz";
//	}

//	public File getTaxCountsFile(String goalname) {
//		return new File(getResultsDir(), goalname + "_counts_" + fastqFile.getName() + ".csv");
//	}

//	public File getKrakenErrFile() {
//		return new File(getResultsDir(), getFastqFile().getName() + getkMserSize() + "_krakenErr.csv");
//	}

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

	public File getKrakenOutDir() {
		return new File(getProjectsDir(), name + "/krakenout");
	}

	public File getTaxIdsFile() {
		return new File(getProjectsDir(), name + "/taxids.txt");
	}

	public File getTaxIdsFilterFile() {
		return new File(getProjectsDir(), name + "/taxidFilter.txt");
	}

	public File getResultsDir() {
		return resFolder;
	}

//	public File getKmerFastqFile() {
//		return getOutputFile(getName(), fastqFile, FileType.FASTQ, true);
//		return new File(getFastqsDir(), getName() + "_k" + getkMserSize() + "_fromFastas.fastq.gz");
//	}

//	public File getKrakenOutFile() {
//		return new File(getKrakenOutDir(), getName() + "_k" + getkMserSize() + "_fromKraken.out.txt.gz");
//	}

//	public File getFilteredKmerFastqFile() {
//		return new File(getFastqsDir(), getName() + "_k" + getkMserSize() + "_fromKraken.fastq.gz");
//	}

//	public File getTrieFile() {
//		return new File(getFiltersDir(), getName() + "_k" + getkMserSize() + "_Trie.ser.gz");
//	}

//	public File getBloomFilterFile() {
//		return new File(getFiltersDir(), getName() + "_k" + getkMserSize() + "_BloomFilter.ser.gz");
//	}
}
