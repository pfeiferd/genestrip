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
import java.util.ArrayList;
import java.util.List;
import java.util.Properties;
import java.util.StringTokenizer;

import org.apache.commons.logging.Log;
import org.metagene.genestrip.genbank.AssemblySummaryReader.FTPEntryQuality;
import org.metagene.genestrip.tax.TaxTree.Rank;
import org.metagene.genestrip.util.GSLogFactory;

public class GSConfig {
	public static final String CONFIG_PROPERTIES = "config.properties";
	public static final String CONFIG_PROPERTIES_2 = "Config.properties";
	public static final String NCBI_FTP_URL = "ftp.ncbi.nih.gov";
	public static final String NCBI_HTTP_BASE_URL = "https://ftp.ncbi.nlm.nih.gov";

	public static final String REF_SEQ_HTTP_BASE_URL = NCBI_HTTP_BASE_URL + "/refseq";
	public static final String REF_SEQ_FTP_URL = NCBI_FTP_URL;

	public static final String TAX_HTTP_BASE_URL = NCBI_HTTP_BASE_URL;
	public static final String TAX_SEQ_FTP_URL = NCBI_FTP_URL;

	public static final String COMPLETE_GENOMES_ONLY = "completeGenomesOnly";
	public static final String IGNORE_MISSING_FASTAS = "ignoreMissingFastas";
	public static final String RANK_COMPLETION_DEPTH = "rankCompletionDepth";
	public static final String MAX_GENOMES_PER_TAXID = "maxGenomesPerTaxid";
	public static final String USE_BLOOM_FILTER_FOR_MATCH = "useBloomFilterForMatch";
	public static final String MAX_READ_TAX_ERROR_COUNT = "maxReadTaxErrorCount";
	public static final String CLASSIFY_READS = "classifyReads";
	public static final String MAX_DUST = "maxDust";
	public static final String SEQ_TYPE = "seqType";
	public static final String REQ_SEQ_LIMIT_FOR_GENBANK = "refSeqLimitForGenbankAccess";
	public static final String MAX_FROM_GENBANK = "maxFromGenBank";
	public static final String FASTA_QUALITIES = "fastaQualities";

	public static enum SeqType {
		GENOMIC, RNA, BOTH;

		public static SeqType byName(String name) {
			for (SeqType r : SeqType.values()) {
				if (r.name().equalsIgnoreCase(name)) {
					return r;
				}
			}
			return null;
		}
	}

	private final File baseDir;
	private final Properties properties;

	public GSConfig(File baseDir) throws IOException {
		this.baseDir = baseDir;
		this.properties = new Properties();
		File configFile = new File(baseDir, CONFIG_PROPERTIES);
		if (!configFile.exists()) {
			configFile = new File(baseDir, CONFIG_PROPERTIES_2);
		}
		Log log = GSLogFactory.getLog("config");
		try {
			if (log.isInfoEnabled()) {
				log.info("Loading config file " + configFile);
			}
			properties.load(new FileInputStream(configFile));
		} catch (IOException e) {
			if (log.isWarnEnabled()) {
				log.warn("Could not read general configuation file '" + configFile + "'. Using defaults.");
			}
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

	public long getMatchLogUpdateCycle() {
		return Long.valueOf(properties.getProperty("matchLogUpdateCycle", "1000000"));
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

	public int getInitialReadSizeBytes() {
		return Integer.valueOf(properties.getProperty("initialReadSizeBytes", "4096"));
	}

	public int getMaxClassificationPaths() {
		return Integer.valueOf(properties.getProperty("maxClassificationPaths", "10"));
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
		return properties.getProperty("krakenBin", "krakenuniq");
	}

	public String getKrakenDB() {
		return properties.getProperty("krakenDB");
	}

	public String getKrakenExecExpr() {
		return properties.getProperty("krakenExecExpr", "{0} -db {1} {2}");
	}

	public boolean isIgnoreMissingFastas() {
		return Boolean.valueOf(properties.getProperty(IGNORE_MISSING_FASTAS, "false"));
	}

	public boolean isClassifyReads() {
		return Boolean.valueOf(properties.getProperty(CLASSIFY_READS, "true"));
	}

	public boolean isUseCompleteGenomesOnly() {
		return Boolean.valueOf(properties.getProperty(COMPLETE_GENOMES_ONLY, "false"));
	}

	public boolean isMatchWithKMerCounts() {
		return Boolean.valueOf(properties.getProperty("matchWithKMerCounts", "false"));
	}

	public int getMaxKMerResCounts() {
		return Integer.valueOf(properties.getProperty("maxKMerResCounts", "200"));
	}

	public double getMaxReadTaxErrorCount() {
		return Double.valueOf(properties.getProperty(MAX_READ_TAX_ERROR_COUNT, "-1"));
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

	public File getFastqDir() {
		return new File(baseDir, "fastq");
	}
	
	public File getFastaDir() {
		return new File(baseDir, "fasta");
	}


	public boolean isUseHttp() {
		return Boolean.valueOf(properties.getProperty("useHttp", "true"));
	}

	public boolean isUseBloomFilterForMatch() {
		return Boolean.valueOf(properties.getProperty(USE_BLOOM_FILTER_FOR_MATCH, "true"));
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

	public Rank getRankCompletionDepth() {
		String s = properties.getProperty(RANK_COMPLETION_DEPTH, null);
		return s == null ? null : Rank.byName(s);
	}

	public SeqType getSeqType() {
		String v = properties.getProperty(GSConfig.SEQ_TYPE, SeqType.GENOMIC.name());
		SeqType dr = SeqType.byName(v);
		if (dr != null) {
			return dr;
		}
		return SeqType.BOTH;
	}

	public int getMaxGenomesPerTaxid() {
		int i = Integer.valueOf(properties.getProperty(MAX_GENOMES_PER_TAXID, "-1"));
		return i <= 0 ? Integer.MAX_VALUE : i;
	}

	public int getRefSeqLimitForGenbankAccess() {
		return Integer.valueOf(properties.getProperty(REQ_SEQ_LIMIT_FOR_GENBANK, "0"));
	}
	
	public int getMaxDust() {
		int i = Integer.valueOf(properties.getProperty(MAX_DUST, "-1"));
		return i <= 0 ? -1 : i;
	}

	public long getNormalizedKMersFactor() {
		return Long.valueOf(properties.getProperty("normalizedKMersFactor", "1000000000"));
	}
	
	public List<FTPEntryQuality> getFastaQualities() {
		String qs = properties.getProperty(FASTA_QUALITIES);
		List<FTPEntryQuality> res = new ArrayList<FTPEntryQuality>();
		if (qs != null) {
			StringTokenizer tokenizer = new StringTokenizer(qs, ",;");
			while (tokenizer.hasMoreTokens()) {
				FTPEntryQuality q = FTPEntryQuality.valueOf(tokenizer.nextToken().trim());
				if (q != null) {
					res.add(q);
				}
			}
		}
		if (res.isEmpty()) {
			res.add(FTPEntryQuality.COMPLETE_LATEST);
			res.add(FTPEntryQuality.CHROMOSOME_LATEST);
		}
		return res;
	}

	public File getGenbankDir() {
		return new File(getCommonDir(), "genbank");		
	}
	
	public int getMaxFromGenbank() {
		String v = properties.getProperty(GSConfig.MAX_FROM_GENBANK, "1");
		if (v != null) {
			return Integer.valueOf(v);
		}
		return -1;
	}	
}
