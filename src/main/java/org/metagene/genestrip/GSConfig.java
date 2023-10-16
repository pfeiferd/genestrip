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
import java.util.Properties;

import org.apache.commons.logging.LogFactory;

public class GSConfig {
	public static final String CONFIG_PROPERTIES = "Config.properties";
	public static final String NCBI_FTP_URL = "ftp.ncbi.nih.gov";
	public static final String NCBI_HTTP_BASE_URL = "https://ftp.ncbi.nlm.nih.gov";

	public static final String REF_SEQ_HTTP_BASE_URL = NCBI_HTTP_BASE_URL + "/refseq";
	public static final String REF_SEQ_FTP_URL = NCBI_FTP_URL;

	public static final String TAX_HTTP_BASE_URL = NCBI_HTTP_BASE_URL;
	public static final String TAX_SEQ_FTP_URL = NCBI_FTP_URL;

	private final File baseDir;
	private final Properties properties;

	public GSConfig(File baseDir) throws IOException {
		this.baseDir = baseDir;
		this.properties = new Properties();
		File configFile = new File(baseDir, CONFIG_PROPERTIES);
		try {
			properties.load(new FileInputStream(configFile));
		} catch (IOException e) {
			LogFactory.getLog("conifg").warn("Could not read configuation file '" + configFile + "'. Using defaults.");
		}
	}

	public String getLogLevel() {
		return properties.getProperty("logLevel", "info");
	}

	public File getBaseDir() {
		return baseDir;
	}

	// Number of threads for goals 'match' and 'filter'.
	// Minimum is 1 (single-threaded), default is -1, which means the number of
	// cores on this computer is used.
	public int getThreads() {
		int threads = Integer.valueOf(properties.getProperty("threads", "-1"));
		if (threads < 0) {
			threads = Runtime.getRuntime().availableProcessors() - 1;
			if (threads < 0) {
				threads = 0;
			}
		}
		return threads;
	}

	public int getThreadQueueSize() {
		return 1000;
	}

	public boolean isCountUniqueKMers() {
		return Boolean.valueOf(properties.getProperty("countUniqueKMers", "true"));
	}

	public boolean isWriteDumpedFastq() {
		return Boolean.valueOf(properties.getProperty("writeDumpedFastq", "false"));
	}

	public boolean isWriteFilteredFastq() {
		return Boolean.valueOf(properties.getProperty("writeFilteredFastq", "false"));
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

	public double getBloomFilterFpp() {
		return Double.valueOf(properties.getProperty("bloomFilterFpp", " 0.00000000001"));
	}

	public String getKrakenBin() {
		return properties.getProperty("krakenBin");
	}

	public String getKrakenDB() {
		return properties.getProperty("krakenDB");
	}

	public String getKrakenExecExpr() {
		return properties.getProperty("krakenExecExpr", "{0} -db {1} {2}");
	}

	public boolean isIgnoreMissingFastas() {
		return Boolean.valueOf(properties.getProperty("ignoreMissingFastas", "false"));
	}

	public boolean isUseCompleteGenomesOnly() {
		return Boolean.valueOf(properties.getProperty("completeGenomesOnly", "false"));
	}

	public boolean isMatchWithKMerCounts() {
		return Boolean.valueOf(properties.getProperty("matchWithKMerCounts", "false"));
	}

	public int getMaxKMerResCounts() {
		return Integer.valueOf(properties.getProperty("maxKMerResCounts", "200"));
	}
	
	public int getMaxReadTaxErrorCount() {
		return Integer.valueOf(properties.getProperty("maxReadTaxErrorCount", "3"));
	}
	
	public String getFtpBaseURL() {
		return properties.getProperty("ftpBaseURL", NCBI_FTP_URL);
	}

	public String getHttpBaseURL() {
		return properties.getProperty("httpBaseURL", NCBI_HTTP_BASE_URL);
	}

	public int getKMerSize() {
		return Integer.valueOf(properties.getProperty("kMerSize", "31"));
	}

	public File getCommonDir() {
		return new File(baseDir, "common");
	}

	public boolean isUseHttp() {
		return Boolean.valueOf(properties.getProperty("useHttp", "true"));
	}

	public File getRefSeqDir() {
		return new File(getCommonDir(), "refseq");
	}

	public String getRefSeqHttpBaseURL() {
		return properties.getProperty("refseqHttpBaseURL", REF_SEQ_HTTP_BASE_URL);
	}

	public String getRefSeqFTPBaseURL() {
		return properties.getProperty("refseqFTPBaseURL", REF_SEQ_FTP_URL);
	}

	public String getTaxHttpBaseURL() {
		return properties.getProperty("taxHttpBaseURL", TAX_HTTP_BASE_URL);
	}

	public String getTaxFTPBaseURL() {
		return properties.getProperty("taxFTPBaseURL", TAX_SEQ_FTP_URL);
	}
}
