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

import org.metagene.genestrip.tax.AssemblySummaryReader.FTPEntryQuality;
import org.metagene.genestrip.tax.TaxTree;
import org.metagene.genestrip.tax.TaxTree.Rank;

public class GSConfig {
	public static final String CONFIG_PROPERTIES = "Config.properties";
	public static final String NCBI_URL = "ftp.ncbi.nih.gov";
	public static final String NCBI_HTTP_BASE_URL = "https://ftp.ncbi.nlm.nih.gov";

	private final File baseDir;
	private final Properties properties;

	public GSConfig(File baseDir) throws IOException {
		this.baseDir = baseDir;
		this.properties = new Properties();
		properties.load(new FileInputStream(new File(baseDir, CONFIG_PROPERTIES)));
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

	public boolean isCountUniqueKmers() {
		return Boolean.valueOf(properties.getProperty("countUniqueKmers", "true"));		
	}
	
	public boolean isIgnoreMissingFastas() {
		return Boolean.valueOf(properties.getProperty("ignoreMissingFastas", "false"));
	}

	public boolean isUseGenBank() {
		return Boolean.valueOf(properties.getProperty("useGenBank", "false"));
	}
	
//	public boolean isUseTrie() {
//		return Boolean.valueOf(properties.getProperty("useTrie", "false"));		
//	}

	public boolean isWriteDumpedFastq() {
		return Boolean.valueOf(properties.getProperty("writeDumpedFastq", "false"));
	}

	public boolean isWriteFilteredFastq() {
		return Boolean.valueOf(properties.getProperty("writeFilteredFastq", "true"));
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

	public double getKMerFastBloomFpp() {
		return Double.valueOf(properties.getProperty("kMerFastBloomFpp", " 0.00000000001"));
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

	public String getSortBin() {
		return properties.getProperty("sortBin", "sort");
	}

	public String getSortExecExpr() {
		return properties.getProperty("sortExecExpr", "{0} -n -t ':' -k 3 {1}");
	}

	public boolean isUseKraken1() {
		return Boolean.valueOf(properties.getProperty("useKraken1", "false"));
	}

	public FTPEntryQuality getFastaQuality() {
		String qs = properties.getProperty("fastaQuality");
		return qs == null ? FTPEntryQuality.COMPLETE_LATEST : FTPEntryQuality.valueOf(qs);
	}

	public String getFtpBaseURL() {
		return properties.getProperty("ftpBaseURL", NCBI_URL);
	}

	public String getHttpBaseURL() {
		return properties.getProperty("httpBaseURL", NCBI_HTTP_BASE_URL);
	}

	public int getKMerSize() {
		return Integer.valueOf(properties.getProperty("kMerSize", "35"));
	}

	public File getCommonDir() {
		return new File(baseDir, "common");
	}

	public File getAdditionalDir() {
		return new File(baseDir, "additional");
	}

	public boolean isUseHttp() {
		return Boolean.valueOf(properties.getProperty("useHttp", "true"));
	}

	public Rank getMaxRankForFilters() {
		return TaxTree.Rank.byName(properties.getProperty("maxRankForFilters", TaxTree.Rank.GENUS.getName()));
	}

	public int getMaxBloomFilterSize() {
		return Integer.valueOf(properties.getProperty("maxBloomFilterSize", Integer.toString(Integer.MAX_VALUE)));
	}
}
