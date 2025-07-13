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

import java.io.PrintStream;
import java.lang.annotation.Annotation;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.StringTokenizer;

import org.metagene.genestrip.genbank.AssemblySummaryReader.AssemblyQuality;
import org.metagene.genestrip.goals.MDDescription;
import org.metagene.genestrip.make.ConfigKey;
import org.metagene.genestrip.make.ConfigParamInfo;
import org.metagene.genestrip.make.ConfigParamInfo.BooleanConfigParamInfo;
import org.metagene.genestrip.make.ConfigParamInfo.DoubleConfigParamInfo;
import org.metagene.genestrip.make.ConfigParamInfo.IntConfigParamInfo;
import org.metagene.genestrip.make.ConfigParamInfo.ListConfigParamInfo;
import org.metagene.genestrip.make.ConfigParamInfo.LongConfigParamInfo;
import org.metagene.genestrip.make.ConfigParamInfo.StringConfigParamInfo;
import org.metagene.genestrip.make.GoalKey;
import org.metagene.genestrip.tax.Rank;

public enum GSConfigKey implements ConfigKey {
	// General
	@MDDescription("Only the log levels `error`, `warn`, `info` and `trace` are used by Genestrip.")
	LOG_LEVEL("logLevel", new LogLevelConfigParamInfo("info")),
	@MDDescription("The number of consumer threads *n* when processing data with respect to the goals `match`, `filter` and also so during the update phase of the `db` goal. "
			+ "There is always one additional thread that reads and uncompresses a corresponding fastq or fasta file (so it is *n + 1* threads in total). "
			+ "When negative, the number of available processors *- 1* is used as *n*. When 0, then the corresponding goals run in single-threaded mode.")
	THREADS("threads", new IntConfigParamInfo(-1, 64, -1), GSGoalKey.DB, GSGoalKey.MATCH, GSGoalKey.MATCHLR, GSGoalKey.FILTER),
	@MDDescription("Whether to show a progress bar on the command line for longer taking process steps.")
	PROGRESS_BAR("progressBar", new BooleanConfigParamInfo(true), GSGoalKey.DB, GSGoalKey.MATCH, GSGoalKey.MATCHLR, GSGoalKey.FILTER),
	@MDDescription("Update period in ms for progress bar (if shown).")
	PROGRESS_BAR_UPDATE("progressBarUpdateMs", new IntConfigParamInfo(100, Integer.MAX_VALUE, 1000), GSGoalKey.DB, GSGoalKey.MATCH, GSGoalKey.MATCHLR, GSGoalKey.FILTER),
	@MDDescription("The number of base pairs *k* for *k*-mers. "
			+ "Changes to this values do *not* affect the memory usage of a database.")
	KMER_SIZE("kMerSize", new IntConfigParamInfo(15, 32, 31), GSGoalKey.DB, GSGoalKey.FILTER, GSGoalKey.MATCH, GSGoalKey.MATCHLR),

	// Genomic data download
	@MDDescription("This base URL will be extended by `/pub/taxonomy/` in order to download the taxonomy file `taxdmp.zip` and by `/genomes/genbank` for files from Genbank.")
	HTTP_BASE_URL("httpBaseURL", new StringConfigParamInfo("https://ftp.ncbi.nlm.nih.gov"), GSGoalKey.DB),
	FTP_BASE_URL("ftpBaseURL", new StringConfigParamInfo("ftp.ncbi.nih.gov"), GSGoalKey.DB),
	@MDDescription("This [mirror](https://www.funet.fi/pub/mirrors/ftp.ncbi.nlm.nih.gov/refseq/) might be considered as an alternative. (No other mirror sites are known.)")
	REF_SEQ_HTTP_BASE_URL("refseq.httpBaseURL", new StringConfigParamInfo("https://ftp.ncbi.nlm.nih.gov/refseq"),
			GSGoalKey.DB),
	REF_SEQ_FTP_BASE_URL("refseq.ftpBaseURL", new StringConfigParamInfo("ftp.ncbi.nih.gov"), GSGoalKey.DB),
	@MDDescription("Use http(s) to download data from NCBI. If `false`, then Genestrip will do anonymous FTP instead (with login and password set to `anonymous`).")
	USE_HTTP("useHttp", new BooleanConfigParamInfo(true), GSGoalKey.DB),
	@MDDescription("If `true`, then a download of files from NCBI will not stop in case a file is missing on the server.")
	IGNORE_MISSING_FASTAS("ignoreMissingFastas", new BooleanConfigParamInfo(false), GSGoalKey.DB),
	@MDDescription("The number of download attempts for a file before giving up.")
	MAX_DOWNLOAD_TRIES("maxDownloadTries", new IntConfigParamInfo(1, 1024, 5), GSGoalKey.DB),
	@MDDescription("Which type of sequence files to include from the RefSeq. RNA files from the RefSeq end with `rna.fna.gz`, whereas genomes end with `genomic.fna.gz`.")
	SEQ_TYPE("seqType", new SeqTypeConfigParamInfo(SeqType.GENOMIC), GSGoalKey.DB),
	@MDDescription("The rank up to which tax ids from `taxids.txt` will be completed by descendants of the taxonomy tree (the set rank included). "
			+ "If not set, the completion will traverse down to the lowest possible levels of the [taxonomy](https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip). "
			+ "Typical values could be `species` or `strain`, but  all values used for assigning ranks in the taxonomy are possible.")
	RANK_COMPLETION_DEPTH("rankCompletionDepth", new RankConfigParamInfo(null), GSGoalKey.DB),
	@MDDescription("If true, then md5 check sums may be skipped by creating and accessing a file named `<file>.md5ok` " +
			"that marks wether the md5 check sum of `<file>` was found to be ok after a previous download of `<file>`.")
	CHECK_SUM_CACHE_FILE("checkSumCacheFile", new BooleanConfigParamInfo(true), GSGoalKey.DB),

	// Limit database size
	@MDDescription("The maximum number of genomes per tax id to be included in the database. "
			+ "Note, that this is an important parameter to control database size, because in some cases, there are thousands of genomic entries per tax id.")
	MAX_GENOMES_PER_TAXID("maxGenomesPerTaxid", new IntConfigParamInfo(1, Integer.MAX_VALUE, Integer.MAX_VALUE),
			GSGoalKey.DB),
	@MDDescription("The limit for the number of *k*-mers per tax id at which adding more *k*-mers for this tax id to the database stops. "
			+ "Note, that this is an important parameter to control database size, because in some cases, there are millions of *k*-mers per tax id.")
	MAX_KMERS_PER_TAXID("maxKMersPerTaxid", new LongConfigParamInfo(0, Long.MAX_VALUE, Long.MAX_VALUE)),
	@MDDescription("The rank for which to consider the parameters `maxGenomesPerTaxid` and `maxKMersPerTaxid`. If `null`, then maximum number of genomes is considered with respect to the direct tax id under which a genome is stored.")
	MAX_GENOMES_PER_TAXID_RANK("maxPerTaxidRank", new RankConfigParamInfo(null)),
	@MDDescription("If `true`, a fastq of fasta file which is downloaded via a URL is always assumed to be g-zipped. Otherwise, it will be considered g-zipped only if the" +
			" file part of the URL ends with `.gz` or `.gzip`.")
	ALWAYS_ASSUME_GZIP("alwaysAssumeGzip", new BooleanConfigParamInfo(true), GSGoalKey.FASTA_MAP, GSGoalKey.FASTQ_MAP),

	// Refseq data selection
	@MDDescription("Whether the [RefSeq](https://ftp.ncbi.nlm.nih.gov/refseq/release/) should be used as the basis for filling the database.")
	REF_SEQ_DB("refseq.filldb", new BooleanConfigParamInfo(true), false, GSGoalKey.FILL_DB),
	@MDDescription("If `true`, then only genomic accessions with the prefixes `AC`, `NC_`, `NZ_` will be considered when filling the database. "
			+ "Otherwise, all genomic accessions will be considered. See [RefSeq accession numbers and molecule types](https://www.ncbi.nlm.nih.gov/books/NBK21091/table/ch18.T.refseq_accession_numbers_and_mole/) for details.")
	COMPLETE_GENOMES_ONLY("refseq.completeGenomesOnly", new BooleanConfigParamInfo(false), GSGoalKey.FILL_DB),
	@MDDescription("Determines whether Genestrip should try to lookup genomic fasta files from Genbank, "
			+ "if the number of corresponding reference genomes from the RefSeq is below the given limit for a requested tax id including its descendants. "
			+ "E.g. `refSeq.limitForGenbankAccess=1` would imply that Genbank is consulted if not a single reference genome is found in the RefSeq for a requested tax id. "
			+ "The default `refSeq.limitForGenbankAccess=0` essentially inactivates this feature."
			+ "In addition, Genbank access is also influenced by the keys `genbank.fastaQualities`, `genbank.maxPerTaxid` and `genbank.referenceOnly` (see below)."
			+ "Note that `refSeq.limitForGenbankAccess` is disregarded if `refseq.filldb=false`.")
	REQ_SEQ_LIMIT_FOR_GENBANK("refSeq.limitForGenbankAccess", new IntConfigParamInfo(0, Integer.MAX_VALUE, 0),
			GSGoalKey.DB),
	@MDDescription("The rank for which to check the limit `refSeq.limitForGenbankAccess`. If `null`, then the limit applies to all requested tax ids and its descendants.")
	REQ_SEQ_LIMIT_FOR_GENBANK_RANK("refSeq.limitForGenbankRank", new RankConfigParamInfo(Rank.SPECIES), GSGoalKey.DB),

	// Genbank data selection
	@MDDescription("Determines the maximum number of fasta files used from Genbank per requested tax id. "
			+ "If this value is <= 0 then all fasta files will be used. "
			+ "Otherwise, if the corresponding number of matching files exceeds `genbank.maxPerTaxid`, then  best ones according to `genbank.fastaQualities` will be retained while adhering to this maximum.")
	MAX_FROM_GENBANK("genbank.maxPerTaxid", new IntConfigParamInfo(-1, Integer.MAX_VALUE, 1), GSGoalKey.DB),
	@MDDescription("Determines the allowed quality levels of fasta files from Genbank. "
			+ "The values must be comma-separated. If a corresponding value is included in the list, "
			+ "then a fasta file for a requested tax id on that quality level will be included, "
			+ "otherwise not (while also respecting the conditions exerted via the keys `refSeq.limitForGenbankAccess` and `genbank.maxPerTaxid`). "
			+ "The quality levels are based on Genbank's [Assembly Summary File](https://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank.txt) (columns `version_status` and `assembly_level`). "
			+ "If the list is empty then no fasta files from Genbank will qualify.")
	FASTA_QUALITIES("genbank.fastaQualities", new FastaQualitiesConfigInfo(Arrays.asList(
			new AssemblyQuality[] { AssemblyQuality.COMPLETE_LATEST, AssemblyQuality.CHROMOSOME_LATEST })), GSGoalKey.DB),
	@MDDescription("Whether only reference genomes are accepted or not. (Reference Genomes must be fetched from GenBank.)")
	REF_GEN_ONLY("genbank.referenceOnly", new BooleanConfigParamInfo(false), false, GSGoalKey.DB),

	// Database generation
	@MDDescription("When generating a database via the goal `db`, any low-complexity *k*-mer with too many "
			+ "repetitive sequences of base pairs may be omitted for storing. To do so, Genestrip employs a simple "
			+ "[genetic dust-filter](https://pubmed.ncbi.nlm.nih.gov/16796549/) for *k*-mers: "
			+ "It assigns a dust value *d* to each *k*-mer, and if *d* >  `maxDust`, then the *k*-mer will not be stored. "
			+ "Given a *k*-mer with *n* repeating base pairs of repeat length *k(1), ... k(n)* with *k(i) > 1*, "
			+ "then *d = fib(k(1)) + ... + fib(k(n))*, where *fib(k(i))* is the Fibonacci number of *k(i)*. "
			+ "E.g., for the *8*-mer `TTTCGGTC`, we have *n = 2* with *k(1) = 3*, *k(2) = 2* and *d = fib(3) + fib(2) = 2 + 1 = 3*. "
			+ "For practical concerns `maxDust = 200` may be suitable. In this case, if *31*-mers were uniformly, randomly generated, "
			+ "then less than 0.0001 % of them would be dropped. If `maxDust = -1`, then dust-filtering is inactive.")
	MAX_DUST("maxDust", new IntConfigParamInfo(-1, Integer.MAX_VALUE, -1), GSGoalKey.DB),
	TEMP_BLOOM_FILTER_FPP("tempBloomFilterFpp", new DoubleConfigParamInfo(0, 1, 0.001d), true, GSGoalKey.TEMPINDEX),
	BLOOM_FILTER_FPP("bloomFilterFpp", new DoubleConfigParamInfo(0, 1, 0.00000000001d), true, GSGoalKey.FILL_DB),
	XOR_BLOOM_HASH("xorBloomHash", new BooleanConfigParamInfo(true)),
	FASTA_LINE_SIZE_BYTES("fastaLineSizeBytes", new IntConfigParamInfo(4096, 65536, 4096), true, GSGoalKey.DB),
	@MDDescription("Perform database update regarding least common ancestors only based on genomes of tax ids as selected for the database generation (and not via all of a super-kingdom's RefSeq genomes).")
	MIN_UPDATE("minUpdate", new BooleanConfigParamInfo(false), false, GSGoalKey.UPDATE_DB),
	@MDDescription("If `true`, then only genomic accessions with the prefixes `AC`, `NC_`, `NZ_` will be considered when updating the database. "
			+ "Otherwise, all genomic accessions will be considered for the update phase. See [RefSeq accession numbers and molecule types](https://www.ncbi.nlm.nih.gov/books/NBK21091/table/ch18.T.refseq_accession_numbers_and_mole/) for details.")
	UPDATE_WITH_COMPLETE_GENOMES_ONLY("refseq.updateWithCompleteGenomesOnly", new BooleanConfigParamInfo(false), GSGoalKey.UPDATE_DB),
	@MDDescription("Wether to delete the temporary database after the final database has been saved or not.")
	REMOVE_TEMP_DB("removeTempDB", new BooleanConfigParamInfo(true), false, GSGoalKey.DB),
	@MDDescription("Stores *k*-mers in steps of `stepSize`. " +
			"E.g. if `stepSize=2` then only every second *k*-mer from a genome is considered for entry into the database.")
	STEP_SIZE("stepSize", new IntConfigParamInfo(1, Integer.MAX_VALUE, 1), GSGoalKey.DB),

	// Match
	@MDDescription("Affects the log level `trace`: Defines after how many reads per fastq file, information on the matching progress is logged. If less than 1, then no progress information is logged.")
	LOG_PROGRESS_UPDATE_CYCLE("logProgressUpdateCycle", new LongConfigParamInfo(0, Long.MAX_VALUE, 1000000), GSGoalKey.MATCH, GSGoalKey.MATCHLR, GSGoalKey.FILTER),
	@MDDescription("Whether to do read classification in the style of Kraken and KrakenUniq. Matching is faster without "
			+ "read classification and the columns `kmers`, `unique kmers` and `max contig length` in resulting CSV files are usually more conclusive anyways - "
			+ "in particular with respect to long reads. When read classification is off, the columns `reads` and `kmers from reads` will be 0 in resulting CSV files.")
	CLASSIFY_READS("classifyReads", new BooleanConfigParamInfo(true), GSGoalKey.MATCH),
	@MDDescription("If `true`, the number of unique *k*-mers will be counted and reported. This requires less than 5% of additional main memory.")
	COUNT_UNIQUE_KMERS("countUniqueKMers", new BooleanConfigParamInfo(true), GSGoalKey.MATCH, GSGoalKey.MATCHLR),
	@MDDescription("If `true`, then the goal `match` writes filtered fastq files in the same way that the goal `filter` does. ")
	WRITE_FILTERED_FASTQ("writeFilteredFastq", new BooleanConfigParamInfo(false), GSGoalKey.MATCH, GSGoalKey.MATCHLR),
	@MDDescription("If `true`, Genestrip will write output files with suffix `.out` in the [Kraken output format](https://ccb.jhu.edu/software/kraken/MANUAL.html#output-format) "
			+ "under `<base dir>/projects/<project_name>/krakenout` covering all reads with at least one matching *k*-mer.")
	WRITE_KRAKEN_STYLE_OUT("writeKrakenStyleOut", new BooleanConfigParamInfo(false), GSGoalKey.MATCH, GSGoalKey.MATCHLR),
	@MDDescription("If `true` a bloom filter will be loaded and used during fastq file analysis (i.e. matching). "
			+ "Using the bloom filter tends to shorten matching time, if the most part of the reads cannot be classified because they contain *no* *k*-mers from the database. "
			+ "Otherwise, using the bloom filter might increase matching time by up to 30%. It also requires more main memory.")
	USE_BLOOM_FILTER_FOR_MATCH("useBloomFilterForMatch", new BooleanConfigParamInfo(true), GSGoalKey.MATCH, GSGoalKey.MATCHLR),
	@MDDescription("The absolute or relative maximum number of *k*-mers that do not need to be in the database for a read to be classified (read error count). "
			+ "If the number is above `maxReadTaxErrorCount`, then the read will not be classified. "
			+ "Otherwise the read will be classified in the same way as [done by Kraken](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2014-15-3-r46/figures/1). "
			+ "If `maxReadTaxErrorCount` is >= 1, then it is interpreted as an absolute number of *k*-mers. "
			+ "Otherwise (and so, if >= 0 and < 1), it is interpreted as the ratio between the *k*-mers not in the database and all *k*-mers of the read. "
			+ "If `maxReadTaxErrorCount` < 0, then the read error count is disregarded, which means that even a single matching *k*-mer will lead to the read's classification.")
	MAX_READ_TAX_ERROR_COUNT("maxReadTaxErrorCount", new DoubleConfigParamInfo(-1, Double.MAX_VALUE, -1),
			GSGoalKey.MATCH, GSGoalKey.MATCHLR),
	@MDDescription("The absolute or relative maximum number of *k*-mers that do not need to be consistent with a read's destined class for the read to be classified (read class error count). "
			+ "If the number is above `maxReadClassErrorCount`, then the read will not be classified. "
			+ "Otherwise the read will be classified in the same way as [done by Kraken](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2014-15-3-r46/figures/1). "
			+ "If `maxReadClassErrorCount` is >= 1, then it is interpreted as an absolute number of *k*-mers. "
			+ "Otherwise (and so, if >= 0 and < 1), it is interpreted as the ratio between the *k*-mers not in the database and all *k*-mers of the read. "
			+ "If `maxReadClassErrorCount` < 0, then the read error count is disregarded, which means that even a single matching *k*-mer will lead to the read's classification.")
	MAX_READ_CLASS_ERROR_COUNT("maxReadClassErrorCount", new DoubleConfigParamInfo(-1, Double.MAX_VALUE, -1),
			GSGoalKey.MATCH, GSGoalKey.MATCHLR),
	@MDDescription("If > 0, the corresponding number of frequencies of the most frequent *k*-mers per tax id will be reported.")
	MAX_KMER_RES_COUNTS("maxKMerResCounts", new IntConfigParamInfo(0, 65536, 0), GSGoalKey.MATCH, GSGoalKey.MATCHLR),
	@MDDescription("If `true`, then an optimized version of the filter / matcher with heavily inlined code is used.")
	USE_INLINED("useInlined", new BooleanConfigParamInfo(false), true, GSGoalKey.MATCH, GSGoalKey.MATCHLR, GSGoalKey.FILTER),
	THREAD_QUEUE_SIZE("threadQueueSize", new IntConfigParamInfo(10, 10000, 1000), true),
	INITIAL_READ_SIZE_BYTES("initialReadSizeBytes", new IntConfigParamInfo(256, 65536, 4096), true),
	MAX_CLASSIFICATION_PATHS("maxClassificationPaths", new IntConfigParamInfo(1, 128, 10), true),

	// Filter
	@MDDescription("If `true`, then `filter` will also generate a fastq file `dumped_...` with all reads not written to the corresponding filtered fastq file.")
	WRITE_DUMPED_FASTQ("writeDumpedFastq", new BooleanConfigParamInfo(false), GSGoalKey.FILTER),
	@MDDescription("The mininum number of a read's *k*-mers to be found in the bloom index such that the read is added to the filtered fastq file. If `minPosCountFilter=0`, then `posRatioFilter` becomes effective.")
	MIN_POS_COUNT_FILTER("minPosCountFilter", new IntConfigParamInfo(0, 1024, 1), GSGoalKey.FILTER),
	@MDDescription("Only effective if `minPosCountFilter=0`: The mininum ratio of a read's *k*-mers to be found in the bloom index such that the read is added to the filtered fastq file.")
	POS_RATIO_FILTER("posRatioFilter", new DoubleConfigParamInfo(0, 1, 0.2), GSGoalKey.FILTER),
	@MDDescription("Whether to process bp probabilities and potentially write them to filtered fastq files. (Takes slightly more memory if `true`.)")
	WITH_PROBS("withProbs", new BooleanConfigParamInfo(false),  GSGoalKey.FILTER, GSGoalKey.MATCH, GSGoalKey.MATCHLR),

	// Kraken
	KRAKEN_BIN("krakenBin", new StringConfigParamInfo("krakenuniq"), true),
	KRAKEN_DB("krakenDB", new StringConfigParamInfo("krakenuniq"), true),
	KRAKEN_EXEC_EXPR("krakenExecExpr", new StringConfigParamInfo("{0} -db {1} {2}"), true);

	private final String name;
	private final ConfigParamInfo<?> param;
	private final boolean internal;
	private final GSGoalKey[] forGoals;

	GSConfigKey(String name, ConfigParamInfo<?> param, GSGoalKey... forGoals) {
		this(name, param, false, forGoals);
	}

	GSConfigKey(String name, ConfigParamInfo<?> param, boolean internal, GSGoalKey... forGoals) {
		this.name = name;
		this.param = param;
		this.internal = internal;
		this.forGoals = forGoals;
	}

	public boolean isInternal() {
		return internal;
	}

	@Override
	public String getName() {
		return name;
	}

	public ConfigParamInfo<?> getInfo() {
		return param;
	}

	public boolean isForGoal(GoalKey forGoal) {
		if (forGoal == null) {
			return true;
		}
		for (GoalKey id : forGoals) {
			if (forGoal.equals(id)) {
				return true;
			}
		}
		return false;
	}
	
	@Override
	public String toString() {
		return getName();
	}

	public static void printMDConfigParamInfo(PrintStream ps, GoalKey filterGoalKey) {
		ps.print('|');
		ps.print("Name");
		ps.print('|');
		ps.print("Type");
		ps.print('|');
		ps.print("Value Range");
		ps.print('|');
		ps.print("Default");
		ps.print('|');
		ps.print("Description");
		ps.print('|');
		ps.print("For Goals");
		ps.print('|');
		ps.println();

		ps.print('|');
		ps.print('-');
		ps.print('|');
		ps.print('-');
		ps.print('|');
		ps.print('-');
		ps.print('|');
		ps.print('-');
		ps.print('|');
		ps.print('-');
		ps.print('|');
		ps.print('-');
		ps.print('|');
		ps.println();

		for (GSConfigKey configKey : GSConfigKey.values()) {
			if (!configKey.isInternal() && configKey.isForGoal(filterGoalKey)) {
				ps.print('|');
				ps.print('`');
				ps.print(configKey.getName());
				ps.print('`');
				ps.print('|');
				ps.print(configKey.getInfo().getTypeDescriptor());
				ps.print('|');
				ps.print(configKey.getInfo().getMDRangeDescriptor());
				ps.print('|');
				ps.print('`');
				ps.print(configKey.getInfo().getMDDefaultValue());
				ps.print('`');
				ps.print('|');
				Annotation[] annotations;
				try {
					annotations = GSConfigKey.class.getField(configKey.name()).getAnnotations();
				} catch (NoSuchFieldException e) {
					throw new RuntimeException(e);
				} catch (SecurityException e) {
					throw new RuntimeException(e);
				}
				for (Annotation annotation : annotations) {
					if (annotation instanceof MDDescription) {
						ps.print(((MDDescription) annotation).value());
						break;
					}
				}
				ps.print('|');
				if (configKey.forGoals.length == 0) {
					ps.print("all");
				} else {
					boolean first = true;
					for (GoalKey key : configKey.forGoals) {
						if (!first) {
							ps.print(", ");
						}
						first = false;
						ps.print('`');
						ps.print(key.getName());
						ps.print('`');
					}
				}
				ps.print('|');
				ps.println();
			}
		}
	}

	public static enum SeqType {
		GENOMIC, RNA, M_RNA, ALL_RNA, ALL;

		public static SeqType byName(String name) {
			for (SeqType r : SeqType.values()) {
				if (r.name().equalsIgnoreCase(name)) {
					return r;
				}
			}
			return null;
		}
	}

	public static class FastaQualitiesConfigInfo extends ListConfigParamInfo<AssemblyQuality> {
		public FastaQualitiesConfigInfo(List<AssemblyQuality> defaults) {
			super(defaults);
		}

		@Override
		protected List<AssemblyQuality> fromString(String qs) {
			List<AssemblyQuality> res = new ArrayList<AssemblyQuality>();
			if (qs != null) {
				StringTokenizer tokenizer = new StringTokenizer(qs, ",;");
				while (tokenizer.hasMoreTokens()) {
					AssemblyQuality q = AssemblyQuality.valueOf(tokenizer.nextToken().trim());
					if (q != null) {
						res.add(q);
					}
				}
			}
			return res;
		}

		@Override
		public String getMDRangeDescriptor() {
			StringBuilder builder = new StringBuilder();
			AssemblyQuality[] types = AssemblyQuality.values();
			for (int i = 0; i < types.length; i++) {
				if (i > 0) {
					builder.append(", ");
				}
				builder.append('`');
				builder.append(types[i].name());
				builder.append('`');
			}
			return builder.toString();
		}

		@Override
		public String getTypeDescriptor() {
			return "list of nominals";
		}

	}

	public static class LogLevelConfigParamInfo extends StringConfigParamInfo {
		private String[] LEVELS = new String[] { "all", "trace", "debug", "info", "warn", "error", "fatal", "off" };

		public LogLevelConfigParamInfo(String defaultValue) {
			super(defaultValue);
		}

		@Override
		protected String fromString(String s) {
			if (s != null) {
				s = s.trim();
				for (int i = 0; i < LEVELS.length; i++) {
					if (LEVELS[i].equalsIgnoreCase(s)) {
						return LEVELS[i];
					}
				}
			}
			return null;
		}

		@Override
		public String getMDRangeDescriptor() {
			StringBuilder builder = new StringBuilder();
			for (int i = 0; i < LEVELS.length; i++) {
				if (i > 0) {
					builder.append(", ");
				}
				builder.append('`');
				builder.append(LEVELS[i]);
				builder.append('`');
			}
			return builder.toString();
		}
	}

	public static class RankConfigParamInfo extends ConfigParamInfo<Rank> {
		public RankConfigParamInfo(Rank defaultValue) {
			super(defaultValue);
		}

		@Override
		public boolean isValueInRange(Object o) {
			return o == null || o instanceof Rank;
		}

		@Override
		protected Rank fromString(String s) {
			return Rank.byName(s);
		}

		@Override
		public String getMDRangeDescriptor() {
			StringBuilder builder = new StringBuilder();
			Rank[] ranks = Rank.values();
			for (int i = 0; i < ranks.length; i++) {
				if (i > 0) {
					builder.append(", ");
				}
				builder.append('`');
				builder.append(ranks[i].getName());
				builder.append('`');
			}
			return builder.toString();
		}

		@Override
		public String getTypeDescriptor() {
			return "nominal";
		}
	}

	public static class SeqTypeConfigParamInfo extends ConfigParamInfo<SeqType> {
		public SeqTypeConfigParamInfo(SeqType defaultValue) {
			super(defaultValue);
		}

		@Override
		protected SeqType fromString(String s) {
			return SeqType.byName(s);
		}

		@Override
		public boolean isValueInRange(Object o) {
			return o instanceof SeqType;
		}

		@Override
		public String getMDRangeDescriptor() {
			StringBuilder builder = new StringBuilder();
			SeqType[] types = SeqType.values();
			for (int i = 0; i < types.length; i++) {
				if (i > 0) {
					builder.append(", ");
				}
				builder.append('`');
				builder.append(types[i].name());
				builder.append('`');
			}
			return builder.toString();
		}

		@Override
		public String getTypeDescriptor() {
			return "nominal";
		}
	}
}
