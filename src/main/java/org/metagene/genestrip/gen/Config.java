package org.metagene.genestrip.gen;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.Properties;

import org.metagene.genestrip.tax.AssemblySummaryReader.FTPEntryQuality;

public class Config {
	public static final String CONFIG_PROPERTIES = "Config.properties";
	public static final String NCBI_URL = "ftp.ncbi.nih.gov";
	public static final String NCBI_HTTP_BASE_URL = "https://ftp.ncbi.nlm.nih.gov";

	private final File baseDir;
	private final Properties properties;

	public Config(File baseDir) throws IOException {
		this.baseDir = baseDir;
		this.properties = new Properties();
		properties.load(new FileInputStream(new File(baseDir, CONFIG_PROPERTIES)));
	}

	public File getBaseDir() {
		return baseDir;
	}
	
	public int getMaxReadSizeBytes() {
		return Integer.valueOf(properties.getProperty("maxReadSizeBytes", "8192"));
	}
	
	public int getMinPosCountFilter() {
		return Integer.valueOf(properties.getProperty("minPosCountFilter", "1"));
	}

	public double getPosRatioFilter() {
		return Double.valueOf(properties.getProperty("posRatioFilter", "0.2"));
	}
	
	public String getKrakenBinFolder() {
		return properties.getProperty("krakenBinFolder");
	}

	public String getKrakenDB() {
		return properties.getProperty("krakenDB");
	}

	public String getKrakenExecExpr() {
		return properties.getProperty("krakenExecExpr", "{0}/kraken2 -db {1} {2}");
	}

	public FTPEntryQuality getFastaQuality() {
		String qs = properties.getProperty("fastaQuality");
		return qs == null ? FTPEntryQuality.COMLPETE_LATEST : FTPEntryQuality.valueOf(qs);
	}

	public String getFtpBaseURL() {
		return properties.getProperty("ftpBaseURL", NCBI_URL);
	}

	public String getHttpBaseURL() {
		return properties.getProperty("httpBaseURL", NCBI_HTTP_BASE_URL);
	}

	public int getkMerSize() {
		return Integer.valueOf(properties.getProperty("kMerSize", "35"));
	}

	public File getCommonBaseDir() {
		return new File(baseDir, "common");
	}

	public boolean isUseHttp() {
		return Boolean.valueOf(properties.getProperty("useHttp", "true"));
	}
}
