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
	private final GSConfig config;
	private final String name;
	private final String krakenDB;
	private final FTPEntryQuality fastaQuality;
	private final int kMserSize;
	private final File fastqFile;
	private final File resFolder;
	private final boolean useKrakenOutFilter;

	public GSProject(GSConfig config, String name, FTPEntryQuality fastaQuality, int kMerSize, String krakenDB,
			File fastqFile, File resFolder, boolean useKrakenOutFilter) {
		this.config = config;
		this.name = name;
		this.fastaQuality = fastaQuality;
		this.kMserSize = kMerSize;
		this.krakenDB = krakenDB;
		this.fastqFile = fastqFile;
		this.resFolder = resFolder != null ? resFolder : new File(getProjectsDir(), name + "/csv");
		this.useKrakenOutFilter = useKrakenOutFilter;
	}

	public File getFastqFile() {
		return fastqFile;
	}

	public File getFastqKrakenOutFile() {
		return new File(getKrakenOutDir(), fastqFile.getName() + getkMserSize() + "_fromKraken.out.txt.gz");
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
		return new File(getFastqsDir(), "Filtered_" + fastqFile.getName() + ".gz");
	}

	public File getDumpFastqFile() {
		return getConfig().isWriteDumpedFastq() ? new File(getResultsDir(), "Dumped_" + fastqFile.getName() + ".gz")
				: null;
	}

	public File getTaxCountsFile(String goalname) {
		return new File(getResultsDir(), goalname + "_Counts_" + fastqFile.getName() + ".csv");
	}

	public File getKrakenErrFile() {
		return new File(getResultsDir(), getFastqFile().getName() + getkMserSize() + "_krakenErr.csv");
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

	public File getKrakenOutDir() {
		return new File(getProjectsDir(), name + "/krakenout");
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
		return new File(getKrakenOutDir(), getName() + "_k" + getkMserSize() + "_fromKraken.out.txt.gz");
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

	public boolean isUseKrakenOutFilter() {
		return useKrakenOutFilter;
	}
}
