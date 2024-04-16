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
import java.io.FileInputStream;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.Properties;

import org.apache.commons.logging.Log;
import org.metagene.genestrip.make.FileDownloadGoal.DownloadProject;
import org.metagene.genestrip.GSConfig.SeqType;
import org.metagene.genestrip.io.StreamingResource;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.tax.TaxTree.Rank;
import org.metagene.genestrip.util.GSLogFactory;

public class GSProject implements DownloadProject {
	public static final String PROJECT_PROPERTIES = "project.properties";
	public static final String PROJECT_PROPERTIES_2 = "Project.properties";

	public enum FileType {
		FASTQ_RES(".fastq"), FASTQ(".fastq"), FASTA(".fasta"), CSV(".csv"), KRAKEN_OUT(".out"), KRAKEN_OUT_RES(".out"),
		SER(".ser"), DB(".kmers.ser"), FILTER(".bloom.ser"), LOG(".log");

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
	private final Properties properties;
	private final String taxids;

	public GSProject(GSConfig config, String name) {
		this(config, name, 31, null, null, null, null, false, null);
	}

	public GSProject(GSConfig config, String name, int kMerSize, String krakenDB, File fastqOrCSVFile, File csvDir,
			File fastqResDir, String taxids) {
		this(config, name, kMerSize, krakenDB, fastqOrCSVFile, csvDir, fastqResDir, true, taxids);
	}

	public GSProject(GSConfig config, String name, int kMerSize, String krakenDB, File fastqOrCSVFile, File csvDir,
			File fastqResDir, boolean loadProps, String taxids) {
		this.config = config;
		this.name = name;
		this.kMserSize = kMerSize;
		this.krakenDB = krakenDB;
		this.fastqOrCSVFile = fastqOrCSVFile;
		this.fastqResDir = fastqResDir;
		this.csvDir = csvDir != null ? csvDir : new File(getProjectDir(), "csv");
		this.taxids = taxids;

		this.properties = new Properties();
		if (loadProps) {
			File configFile = new File(getProjectDir(), PROJECT_PROPERTIES);
			if (!configFile.exists()) {
				configFile = new File(getProjectDir(), PROJECT_PROPERTIES_2);
			}
			Log log = GSLogFactory.getLog("project");
			try {
				if (log.isInfoEnabled()) {
					log.info("Loading project file " + configFile);
				}
				properties.load(new FileInputStream(configFile));
			} catch (IOException e) {
				if (log.isWarnEnabled()) {
					log.warn("Could not read project configuation file '" + configFile + "'. Using defaults.");
				}
			}
		}
	}
	
	public String getTaxids() {
		return taxids;
	}

	public File getFastqOrCSVFile() {
		return fastqOrCSVFile;
	}

	public File getDirForType(FileType type) {
		switch (type) {
		case FASTQ_RES:
			return getFastqResDir();
		case FASTQ:
			return getFastqDir();
		case FASTA:
			return getFastaDir();
		case CSV:
			return getResultsDir();
		case KRAKEN_OUT:
		case KRAKEN_OUT_RES:
			return getKrakenOutDir();
		case SER:
		case DB:
		case FILTER:
			return getDBDir();
		case LOG:
			return getLogDir();
		default:
			throw new IllegalArgumentException("Illegal FileType: " + type);
		}
	}

	public File getOutputFile(String goal, FileType type) {
		return getOutputFile(goal, type, true);
	}

	public File getOutputFile(String goal, FileType type, boolean gzip) {
		return getOutputFile(goal, null, null, type, gzip);
	}

	public File getOutputFile(String goal, String baseFile, FileType type) {
		return getOutputFile(goal, null, baseFile, type, true);
	}

	public File getOutputFile(String goal, String key, String baseFile, FileType type, boolean gzip) {
		String baseName = baseFile == null ? "" : baseFile;
		if (baseName.endsWith(".gz")) {
			baseName = baseName.substring(0, baseName.length() - 3);
		} else if (baseName.endsWith(".gzip")) {
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

		return new File(getDirForType(type), getOutputFilePrefix(goal) + getOutputFileGoalPrefix(goal, key) + baseName
				+ type.getSuffix() + (gzip ? ".gz" : ""));
	}

	protected String getOutputFilePrefix(String goal) {
		return getName() + "_";
	}

	protected String getOutputFileGoalPrefix(String goal, String key) {
		if (goal == null) {
			return key == null ? "" : key;
		} else {
			return key == null ? goal : goal + "_" + key;
		}
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

	public String getName() {
		return name;
	}

	public File getProjectsDir() {
		return new File(config.getBaseDir(), "projects");
	}

	public File getProjectDir() {
		return new File(getProjectsDir(), name);
	}

	public File getFastaDir() {
		return new File(getProjectDir(), "fasta");
	}

	public File getFastqDir() {
		return new File(getProjectDir(), "fastq");
	}

	public File getFastqResDir() {
		return fastqResDir == null ? getFastqDir() : fastqResDir;
	}

	public File getDBDir() {
		return new File(getProjectDir(), "db");
	}

	public File getKrakenOutDir() {
		return new File(getProjectDir(), "krakenout");
	}

	public File getTaxIdsFile() {
		return new File(getProjectDir(), "taxids.txt");
	}

	public File getAddtionalFile() {
		return new File(getProjectDir(), "additional.txt");
	}

	public File getCategoriesFile() {
		return new File(getProjectDir(), "categories.txt");
	}

	public File getResultsDir() {
		return csvDir;
	}

	public File getLogDir() {
		return new File(getProjectDir(), "log");
	}

	public boolean isUseBloomFilterForMatch() {
		String v = properties.getProperty(GSConfig.USE_BLOOM_FILTER_FOR_MATCH);
		if (v != null) {
			return Boolean.valueOf(v);

		}
		return getConfig().isUseBloomFilterForMatch();
	}

	public boolean isUseCompleteGenomesOnly() {
		String v = properties.getProperty(GSConfig.COMPLETE_GENOMES_ONLY);
		if (v != null) {
			return Boolean.valueOf(v);

		}
		return getConfig().isUseCompleteGenomesOnly();
	}

	@Override
	public boolean isIgnoreMissingFiles() {
		String v = properties.getProperty(GSConfig.IGNORE_MISSING_FASTAS);
		if (v != null) {
			return Boolean.valueOf(v);

		}
		return getConfig().isIgnoreMissingFastas();
	}

	public boolean isClassifyReads() {
		String v = properties.getProperty(GSConfig.CLASSIFY_READS);
		if (v != null) {
			return Boolean.valueOf(v);

		}
		return getConfig().isClassifyReads();
	}

	public Rank getRankCompletionDepth() {
		String v = properties.getProperty(GSConfig.RANK_COMPLETION_DEPTH);
		if (v != null) {
			return Rank.byName(v);
		}
		return getConfig().getRankCompletionDepth();
	}

	public SeqType getSeqType() {
		String v = properties.getProperty(GSConfig.SEQ_TYPE);
		if (v != null) {
			SeqType dr = SeqType.byName(v);
			if (dr != null) {
				return dr;
			}
		}
		return getConfig().getSeqType();
	}
	
	public int getMaxGenomesPerTaxid() {
		String v = properties.getProperty(GSConfig.MAX_GENOMES_PER_TAXID);
		if (v != null) {
			int i = Integer.valueOf(v);
			return i <= 0 ? Integer.MAX_VALUE : i;
		}
		return getConfig().getMaxGenomesPerTaxid();
	}

	public double getMaxReadTaxErrorCount() {
		String v = properties.getProperty(GSConfig.MAX_READ_TAX_ERROR_COUNT);
		if (v != null) {
			return Double.valueOf(v);
		}
		return getConfig().getMaxReadTaxErrorCount();
	}

	public int getMaxDust() {
		String v = properties.getProperty(GSConfig.MAX_DUST);
		if (v != null) {
			int i = Integer.valueOf(v);
			return i <= 0 ? -1 : i;
		}
		return getConfig().getMaxDust();
	}

	public File getFilterFile(Goal<GSProject> goal) {
		return getOutputFile(goal.getName(), FileType.FILTER);
	}

	public Map<String, List<StreamingResource>> getKeyToFastqs() {
		return null;
	}
}
