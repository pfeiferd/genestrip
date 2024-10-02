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
import java.io.InputStream;
import java.io.UnsupportedEncodingException;
import java.net.URLEncoder;
import java.nio.file.FileSystems;
import java.nio.file.Path;
import java.nio.file.PathMatcher;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Properties;

import org.metagene.genestrip.io.StreamingResource;
import org.metagene.genestrip.make.ConfigKey;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.Project;
import org.metagene.genestrip.util.GSLogFactory;

public class GSProject extends Project {
	public static final String CONFIG_PROPERTIES = "config.properties";
	public static final String CONFIG_PROPERTIES_2 = "Config.properties";

	public enum FileType {
		FASTQ_RES(".fastq"), FASTQ(".fastq"), FASTA(".fasta"), CSV(".csv"), KRAKEN_OUT(".out"), KRAKEN_OUT_RES(".out"),
		SER(".ser"), DB(".zip"), FILTER(".ser"), LOG(".log");

		private final String suffix;

		private FileType(String suffix) {
			this.suffix = suffix;
		}

		public String getSuffix() {
			return suffix;
		}
	}

	private final GSCommon common;
	private final String key;
	private final String[] fastqResources;
	private final String fastqMapFile;
	private final File csvDir;
	private final File fastqResDir;
	private final Properties[] properties;
	private final String taxids;
	private final String dbPath;

	private boolean downloadFastqs;
	private boolean downloadFastqsToCommon;

	public GSProject(GSCommon config, String name) {
		this(config, name, null, null, null, null, null, true, null, null, null, null, false);
	}

	public GSProject(GSCommon config, String name, boolean quiet) {
		this(config, name, null, null, null, null, null, true, null, null, null, null, quiet);
	}

	public GSProject(GSCommon config, String name, String key, String[] fastqFiles) {
		this(config, name, key, fastqFiles, null, null, null, true, null, null, null, null, false);
	}

	public GSProject(GSCommon config, String name, String key, String[] fastqFiles, boolean quiet) {
		this(config, name, key, fastqFiles, null, null, null, true, null, null, null, null, quiet);
	}

	public GSProject(GSCommon config, String name, String csvFile) {
		this(config, name, null, null, csvFile, null, null, true, null, null, null, null, false);
	}

	public GSProject(GSCommon config, String name, String key, String[] fastqFiles, String csvFile, File csvDir,
			File fastqResDir, boolean loadProps, String taxids, Properties commandLineProps, GSGoalKey forGoal,
			String dbPath, boolean quietInit) {
		super(name);
		this.common = config;
		this.fastqResDir = fastqResDir;
		this.key = key;
		this.fastqResources = fastqFiles;
		this.fastqMapFile = csvFile;
		this.csvDir = csvDir != null ? csvDir : new File(getProjectDir(), "csv");
		this.taxids = taxids;
		this.dbPath = dbPath;

		if (commandLineProps != null) {
			initConfigParams(commandLineProps); // Gotta do this here already, e.g. to impact logging early on.
			configureLogger();
		}

		properties = new Properties[3];
		properties[0] = commandLineProps == null ? new Properties() : commandLineProps;
		properties[1] = loadConfigProperties(getProjectDir(), quietInit);
		properties[2] = loadConfigProperties(getCommon().getBaseDir(), quietInit);

		initConfigParams(properties);
		configureLogger();

		for (Properties props : properties) {
			checkConfigProperties(props, forGoal);
		}
		if (!quietInit) {
			logParamMap();
		}
	}

	protected Properties loadConfigProperties(File dir, boolean quiet) {
		Properties properties = new Properties();
		File configFile = new File(dir, CONFIG_PROPERTIES);
		if (!configFile.exists()) {
			configFile = new File(dir, CONFIG_PROPERTIES_2);
		}
		if (!quiet && getLogger().isInfoEnabled()) {
			getLogger().info("Loading config file '" + configFile + "'.");
		}
		try (InputStream is = new FileInputStream(configFile)) {
			properties.load(is);
		} catch (IOException e) {
			if (!quiet && getLogger().isWarnEnabled()) {
				getLogger().warn("Could not read config file '" + configFile + "'.");
			}
		}
		return properties;
	}

	protected void configureLogger() {
		GSLogFactory.getInstance().setLogLevel(stringConfigValue(GSConfigKey.LOG_LEVEL));
	}

	@Override
	protected ConfigKey[] getConfigKeys() {
		return GSConfigKey.values();
	}

	public String getKey() {
		return key;
	}

	public String getDBPath() {
		return dbPath;
	}

	public String[] getFastqResources() {
		return fastqResources;
	}
	
	// For override...
	public List<StreamingResource> getStreamingFastqResources() {
		return null;
	}

	public String getFastqMapFile() {
		return fastqMapFile;
	}

	public String getTaxids() {
		return taxids;
	}

	public boolean isDownloadFastqs() {
		return downloadFastqs;
	}

	public boolean isDownloadFastqsToCommon() {
		return downloadFastqsToCommon;
	}

	public void setDownloadFastqs(boolean downloadFastqs) {
		this.downloadFastqs = downloadFastqs;
	}

	public void setDownloadFastqsToCommon(boolean downloadFastqsToCommon) {
		this.downloadFastqsToCommon = downloadFastqsToCommon;
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
		return getOutputFile(getDirForType(type), goal, key, baseFile, type, gzip);
	}

	public File getOutputFile(File dir, String goal, String key, String baseFile, FileType type, boolean gzip) {
		return getOutputFile(dir, getOutputFilePrefix(goal), goal, key, baseFile, type, gzip);
	}

	public File getOutputFile(File dir, String project, String goal, String key, String baseFile, FileType type,
			boolean gzip) {
		String baseName = baseFile == null ? "" : getFileBaseName(baseFile);
		if (baseName.startsWith(getName() + "_")) {
			baseName = baseName.substring(getName().length() + 1);
		}
		if (!baseName.isEmpty()) {
			baseName = "_" + baseName;
		}

		return new File(dir, (project == null ? "" : project) + getOutputFileGoalPrefix(goal, key) + baseName
				+ type.getSuffix() + (gzip ? ".gz" : ""));
	}

	public String getFileBaseName(String projectFileName) {
		String baseName = projectFileName;
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
		return baseName;
	}

	protected String getOutputFilePrefix(String goal) {
		return getName() + "_";
	}

	protected String getOutputFileGoalPrefix(String goal, String key) {
		if (key != null) {
			try {
				key = URLEncoder.encode(key, "UTF-8");
				if (key.length() > 256) {
					key = key.substring(0, 256);
				}
			} catch (UnsupportedEncodingException e) {
				throw new RuntimeException(e);
			}
		}
		if (goal == null) {
			return key == null ? "" : key;
		} else {
			return key == null ? goal : goal + "_" + key;
		}
	}

	public GSCommon getCommon() {
		return common;
	}

	public File getProjectsDir() {
		return new File(getCommon().getBaseDir(), "projects");
	}

	public File getProjectDir() {
		return new File(getProjectsDir(), getName());
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

	public File getDBFile() {
		return getOutputFile(GSGoalKey.DB.getName(), FileType.DB, false);
	}

	public File getDBInfoFile() {
		return getOutputFile(GSGoalKey.DBINFO.getName(), FileType.CSV, false);
	}

	public File getKrakenOutDir() {
		return new File(getProjectDir(), "krakenout");
	}

	public File getTaxIdsFile() {
		return new File(getProjectDir(), "taxids.txt");
	}

	public File getAdditionalFile() {
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

	public File getFilterFile(Goal<GSProject> goal) {
		return getOutputFile(goal.getKey().getName(), FileType.FILTER);
	}

	public File fastaFileFromPath(String fastaFilePath) {
		if (fastaFilePath != null) {
			File fasta = new File(fastaFilePath);
			if (fasta.exists()) {
				return fasta;
			}
			fasta = new File(getFastaDir(), fastaFilePath);
			if (fasta.exists()) {
				return fasta;
			}
			fasta = new File(getCommon().getFastaDir(), fastaFilePath);
			if (fasta.exists()) {
				return fasta;
			}
		}
		return null;
	}

	public File fastqFileFromPath(String fastqFilePath) {
		if (fastqFilePath != null) {
			File fastq = new File(fastqFilePath);
			if (fastq.exists()) {
				return fastq;
			}
			fastq = new File(getFastqDir(), fastqFilePath);
			if (fastq.exists()) {
				return fastq;
			}
			fastq = new File(getCommon().getFastqDir(), fastqFilePath);
			if (fastq.exists()) {
				return fastq;
			}
		}
		return null;
	}

	public List<File> fastqFilesFromPath(String fastqFilePath) {
		if (fastqFilePath != null) {
			File fastq = new File(fastqFilePath);
			if (fastq.exists()) {
				return Collections.singletonList(fastq);
			}
			fastq = new File(getFastqDir(), fastqFilePath);
			if (fastq.exists()) {
				return Collections.singletonList(fastq);
			}
			fastq = new File(getCommon().getFastqDir(), fastqFilePath);
			if (fastq.exists()) {
				return Collections.singletonList(fastq);
			}
			List<File> list = null;
			list = findFilesByGlobPattern(null, fastqFilePath);
			if (list != null) {
				return list;
			}
			list = findFilesByGlobPattern(getFastqDir(), fastqFilePath);
			if (list != null) {
				return list;
			}
			list = findFilesByGlobPattern(getCommon().getFastqDir(), fastqFilePath);
			if (list != null) {
				return list;
			}
		}
		return null;
	}

	protected List<File> findFilesByGlobPattern(File rootDir, String pattern) {
		List<File> res = null;
		File dir;
		if (rootDir == null) {
			dir = new File(pattern).getParentFile();
		}
		else {
			dir = new File(rootDir, pattern).getParentFile();
		}
		if (dir != null) {
			File[] files = dir.listFiles();
			if (files != null && files.length > 0) {
				pattern = new File(pattern).getName();
				PathMatcher matcher = FileSystems.getDefault().getPathMatcher("glob:" + pattern);
				for (File file : files) {
					Path path = file.toPath();
					if (matcher.matches(rootDir == null ? path : path.getFileName())) {
						if (res == null) {
							res = new ArrayList<File>();
						}
						res.add(file);
					}
				}
			}
		}
		return res;
	}
}
