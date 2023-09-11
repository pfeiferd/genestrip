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
import java.io.Serializable;
import java.util.List;

import org.metagene.genestrip.make.FileDownloadGoal.DownloadProject;
import org.metagene.genestrip.store.KMerSortedArray;
import org.metagene.genestrip.store.KMerStore;
import org.metagene.genestrip.store.KMerStoreFactory;
import org.metagene.genestrip.tax.AssemblySummaryReader.FTPEntryQuality;

public class GSProject implements DownloadProject,  KMerStoreFactory {
	public enum FileType {
		FASTQ_RES(".fastq"), FASTQ(".fastq"), FASTA(".fasta"), CSV(".csv"), KRAKEN_OUT(".out"), KRAKEN_OUT_RES(".out"), SER(".ser");

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
	private final int kMserSize;
	private final File fastqOrCSVFile;
	private final File csvDir;
	private final File fastqResDir;

	public GSProject(GSConfig config, String name, int kMerSize, String krakenDB,
			File fastqOrCSVFile, File csvDir, File fastqResDir) {
		this.config = config;
		this.name = name;
		this.kMserSize = kMerSize;
		this.krakenDB = krakenDB;
		this.fastqOrCSVFile = fastqOrCSVFile;
		this.fastqResDir = fastqResDir;
		this.csvDir = csvDir != null ? csvDir : new File(getProjectsDir(), name + "/csv");
	}

	public File getFastqOrCSVFile() {
		return fastqOrCSVFile;
	}
	
	public File getDirForType(FileType type) {
		switch (type) {
		case FASTQ_RES:
			return getFastqsResDir();
		case FASTQ:
			return getFastqsDir();
		case FASTA:
			return getFastasDir();
		case CSV:
			return getResultsDir();
		case KRAKEN_OUT:
		case KRAKEN_OUT_RES:
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
			String suffix = ft.getSuffix();
			if (baseName.endsWith(suffix)) {
				baseName = baseName.substring(0, baseName.length() - suffix.length());
			}			
		}		
		if (baseName.startsWith(getName() + "_")) {
			baseName = baseName.substring(getName().length() + 1);
		}
		if (!baseName.isEmpty()) {
			baseName = "_" + baseName;
		}
		
		return new File(getDirForType(type), getName() + "_" + goal + baseName + type.getSuffix() + (gzip ? ".gz" : ""));
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

	public GSConfig getConfig() {
		return config;
	}

	public String getKrakenDB() {
		return krakenDB == null ? config.getKrakenDB() : krakenDB;
	}

	public File getBaseDir() {
		return config.getBaseDir();
	}

	public int getKMserSize() {
		return kMserSize <= 0 ? config.getKMerSize() : kMserSize;
	}

	public List<FTPEntryQuality> getFastaQuality() {
		return config.getFastaQualities();
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
	
	public File getFastqsResDir() {
		return fastqResDir == null ? getFastqsDir() : fastqResDir;
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

	public File getTaxIdsFilterFile2() {
		return new File(getProjectsDir(), name + "/taxidFilter2.txt");
	}
	
	public File getResultsDir() {
		return csvDir;
	}
	
	public KMerStoreFactory getKMerStoreFactory() {
		return this;
	}
	
	@Override
	public <V extends Serializable> KMerStore<V> createKMerStore(Class<V> clazz, Object... params) {
//		if (config.isUseTrie()) {
//			return new KMerTrie<V>(2, getKMserSize(), false);
//		}
//		else {
			return new KMerSortedArray<V>(getKMserSize(), 0.000000001, null, false, true);
//		}
	}
}
