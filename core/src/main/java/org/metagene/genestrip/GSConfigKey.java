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
import java.util.*;

import org.metagene.genestrip.bloom.BlockedKMerBloomFilter;
import org.metagene.genestrip.store.RadixKMerStore;
import org.metagene.genestrip.genbank.AssemblySummaryReader.AssemblyQuality;
import org.metagene.genestrip.make.*;
import org.metagene.genestrip.make.ConfigParamInfo.BooleanConfigParamInfo;
import org.metagene.genestrip.make.ConfigParamInfo.DoubleConfigParamInfo;
import org.metagene.genestrip.make.ConfigParamInfo.IntConfigParamInfo;
import org.metagene.genestrip.make.ConfigParamInfo.ListConfigParamInfo;
import org.metagene.genestrip.make.ConfigParamInfo.LongConfigParamInfo;
import org.metagene.genestrip.make.ConfigParamInfo.StringConfigParamInfo;
import org.metagene.genestrip.tax.Rank;

/**
 * Enumeration of all Genestrip configuration parameters. Each constant binds a configuration name
 * to its {@link ConfigParamInfo} (type, value range and default) and to the goals it applies to.
 * Constants annotated with {@link MDDescription} contribute to the generated documentation.
 */
public enum GSConfigKey implements ConfigKey {
	// General
	/** The log level used by Genestrip. */
	@MDDescription("Only the log levels `error`, `warn`, `info` and `trace` are used by Genestrip.")
	LOG_LEVEL("logLevel", new LogLevelConfigParamInfo("info")),
	/** Number of consumer threads used for matching, filtering and the database update phase. */
	@MDDescription("The number of consumer threads *n* when processing data with respect to the goals `match`, `filter` and also so during the update phase of the `db` goal. "
			+ "There is always one additional thread that reads and uncompresses a corresponding fastq or fasta file (so it is *n + 1* threads in total). "
			+ "When negative, the number of available processors *- 1* is used as *n*. When 0, then the corresponding goals run in single-threaded mode.")
	THREADS("threads", new IntConfigParamInfo(-1, 64, -1), GSGoalKey.DB, GSGoalKey.MATCH, GSGoalKey.MATCHLR, GSGoalKey.FILTER),
	/** Whether to show a progress bar for longer running steps. */
	@MDDescription("Whether to show a progress bar on the command line for longer taking process steps.")
	PROGRESS_BAR("progressBar", new BooleanConfigParamInfo(true), GSGoalKey.DB, GSGoalKey.MATCH, GSGoalKey.MATCHLR, GSGoalKey.FILTER),
	/** Update period in milliseconds for the progress bar. */
	@MDDescription("Update period in ms for progress bar (if shown).")
	PROGRESS_BAR_UPDATE("progressBarUpdateMs", new IntConfigParamInfo(100, Integer.MAX_VALUE, 1000), GSGoalKey.DB, GSGoalKey.MATCH, GSGoalKey.MATCHLR, GSGoalKey.FILTER),
	/** The k-mer length k in base pairs. */
	@MDDescription("The number of base pairs *k* for *k*-mers. "
			+ "Changes to this values do *not* affect the memory usage of a database.")
	// Max is 31, not 32: at k=32 a k-mer fills all 64 bits, so the all-T k-mer equals the -1L
	// "invalid k-mer" sentinel used throughout CGAT/the matcher (causing a hang) and longToKMer
	// would see negative values. k<=31 keeps the top bits free for the sentinel.
	KMER_SIZE("kMerSize", new IntConfigParamInfo(15, 31, 31), GSGoalKey.DB, GSGoalKey.FILTER, GSGoalKey.MATCH, GSGoalKey.MATCHLR),
	/** Extract key that read descriptors must match for extraction. */
	@MDDescription("Extract key for read descriptors. The beginning of a descriptor must match this key after the '@' for the read to be written.")
	EXTRACT_KEY("extractKey", new StringConfigParamInfo(""), GSGoalKey.EXTRACT),

	// Genomic data download
	/** Base URL for downloading the taxonomy and Genbank files. */
	@MDDescription("This base URL will be extended by `/pub/taxonomy/` in order to download the taxonomy file `taxdmp.zip` and by `/genomes/genbank` for files from Genbank.")
	HTTP_BASE_URL("httpBaseURL", new StringConfigParamInfo("https://ftp.ncbi.nlm.nih.gov"), GSGoalKey.DB),
	/** FTP base URL for NCBI downloads. */
	FTP_BASE_URL("ftpBaseURL", new StringConfigParamInfo("ftp.ncbi.nih.gov"), GSGoalKey.DB),
	/** HTTP base URL for downloading RefSeq data. */
	@MDDescription("This [mirror](https://www.funet.fi/pub/mirrors/ftp.ncbi.nlm.nih.gov/refseq/) might be considered as an alternative. (No other mirror sites are known.)")
	REF_SEQ_HTTP_BASE_URL("refseq.httpBaseURL", new StringConfigParamInfo("https://ftp.ncbi.nlm.nih.gov/refseq"),
			GSGoalKey.DB),
	/** FTP base URL for downloading RefSeq data. */
	REF_SEQ_FTP_BASE_URL("refseq.ftpBaseURL", new StringConfigParamInfo("ftp.ncbi.nih.gov"), GSGoalKey.DB),
	/** Whether to download NCBI data over HTTP(S) rather than FTP. */
	@MDDescription("Use http(s) to download data from NCBI. If `false`, then Genestrip will do anonymous FTP instead (with login and password set to `anonymous`).")
	USE_HTTP("useHttp", new BooleanConfigParamInfo(true), GSGoalKey.DB),
	/** Whether to continue downloading when a file is missing on the server. */
	@MDDescription("If `true`, then a download of files from NCBI will not stop in case a file is missing on the server.")
	IGNORE_MISSING_FASTAS("ignoreMissingFastas", new BooleanConfigParamInfo(false), GSGoalKey.DB),
	/** Number of download attempts per file before giving up. */
	@MDDescription("The number of download attempts for a file before giving up.")
	MAX_DOWNLOAD_TRIES("maxDownloadTries", new IntConfigParamInfo(1, 1024, 5), GSGoalKey.DB),
	/** Which type of RefSeq sequence files to include. */
	@MDDescription("Which type of sequence files to include from the RefSeq. RNA files from the RefSeq end with `rna.fna.gz`, whereas genomes end with `genomic.fna.gz`.")
	SEQ_TYPE("seqType", new SeqTypeConfigParamInfo(SeqType.GENOMIC), GSGoalKey.DB),
	/** Rank up to which requested tax ids are completed by taxonomy descendants. */
	@MDDescription("The rank up to which tax ids from `taxids.txt` will be completed by descendants of the taxonomy tree (the set rank included). "
			+ "If not set, the completion will traverse down to the lowest possible levels of the [taxonomy](https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip). "
			+ "Typical values could be `species` or `strain`, but  all values used for assigning ranks in the taxonomy are possible.")
	RANK_COMPLETION_DEPTH("rankCompletionDepth", new RankConfigParamInfo(null), GSGoalKey.DB),
	/** Whether md5 checksums may be skipped via a cached {@code .md5ok} marker file. */
	@MDDescription("If true, then md5 check sums may be skipped by creating and accessing a file named `<file>.md5ok` " +
			"that marks whether the md5 check sum of `<file>` was found to be ok after a previous download of `<file>`.")
	CHECK_SUM_CACHE_FILE(FileDownloadGoal.CHECK_SUM_CACHE_FILE.getName(), FileDownloadGoal.CHECK_SUM_CACHE_FILE.getInfo(), GSGoalKey.DB),

	// Limit database size
	/** Maximum number of genomes per tax id included in the database. */
	@MDDescription("The maximum number of genomes per tax id to be included in the database. "
			+ "Note, that this is an important parameter to control database size, because in some cases, there are thousands of genomic entries per tax id.")
	MAX_GENOMES_PER_TAXID("maxGenomesPerTaxid", new IntConfigParamInfo(1, Integer.MAX_VALUE, Integer.MAX_VALUE),
			GSGoalKey.DB),
	/** Maximum number of k-mers stored per tax id. */
	@MDDescription("The limit for the number of *k*-mers per tax id at which adding more *k*-mers for this tax id to the database stops. "
			+ "Note, that this is an important parameter to control database size, because in some cases, there are millions of *k*-mers per tax id.")
	MAX_KMERS_PER_TAXID("maxKMersPerTaxid", new LongConfigParamInfo(0, Long.MAX_VALUE, Long.MAX_VALUE)),
	/** Rank at which the per-tax-id genome and k-mer limits apply. */
	@MDDescription("The rank for which to consider the parameters `maxGenomesPerTaxid` and `maxKMersPerTaxid`. If `null`, then maximum number of genomes is considered with respect to the direct tax id under which a genome is stored.")
	MAX_GENOMES_PER_TAXID_RANK("maxPerTaxidRank", new RankConfigParamInfo(null)),
	/** Whether downloaded fastq/fasta files are always assumed to be gzipped. */
	@MDDescription("If `true`, a fastq or fasta file which is downloaded via a URL is always assumed to be g-zipped. Otherwise, it will be considered g-zipped only if the" +
			" file part of the URL ends with `.gz` or `.gzip`.")
	ALWAYS_ASSUME_GZIP("alwaysAssumeGzip", new BooleanConfigParamInfo(true), GSGoalKey.FASTA_MAP, GSGoalKey.FASTQ_MAP),

	// Refseq data selection
	/** Whether the RefSeq is used as the basis for filling the database. */
	@MDDescription("Whether the [RefSeq](https://ftp.ncbi.nlm.nih.gov/refseq/release/) should be used as the basis for filling the database.")
	REF_SEQ_DB("refseq.filldb", new BooleanConfigParamInfo(true), false, GSGoalKey.FILL_DB),
	/** Whether to consider only complete-genome accession prefixes when filling the database. */
	@MDDescription("If `true`, then only genomic accessions with the prefixes `AC`, `NC_`, `NZ_` will be considered when filling the database. "
			+ "Otherwise, all genomic accessions will be considered. See [RefSeq accession numbers and molecule types](https://www.ncbi.nlm.nih.gov/books/NBK21091/table/ch18.T.refseq_accession_numbers_and_mole/) for details.")
	COMPLETE_GENOMES_ONLY("refseq.completeGenomesOnly", new BooleanConfigParamInfo(false), GSGoalKey.FILL_DB),
	/** Threshold below which Genbank is consulted for additional genomes. */
	@MDDescription("Determines whether Genestrip should try to lookup genomic fasta files from Genbank, "
			+ "if the number of corresponding reference genomes from the RefSeq is below the given limit for a requested tax id including its descendants. "
			+ "E.g. `refSeq.limitForGenbankAccess=1` would imply that Genbank is consulted if not a single reference genome is found in the RefSeq for a requested tax id. "
			+ "The default `refSeq.limitForGenbankAccess=0` essentially inactivates this feature. "
			+ "In addition, Genbank access is also influenced by the keys `genbank.fastaQualities`, `genbank.maxPerTaxid` and `genbank.referenceOnly` (see below). "
			+ "Note that `refSeq.limitForGenbankAccess` is disregarded if `refseq.filldb=false`.")
	REQ_SEQ_LIMIT_FOR_GENBANK("refSeq.limitForGenbankAccess", new IntConfigParamInfo(0, Integer.MAX_VALUE, 0),
			GSGoalKey.DB),
	/** Rank at which the Genbank access limit is checked. */
	@MDDescription("The rank for which to check the limit `refSeq.limitForGenbankAccess`. If `null`, then the limit applies to all requested tax ids and its descendants.")
	REQ_SEQ_LIMIT_FOR_GENBANK_RANK("refSeq.limitForGenbankRank", new RankConfigParamInfo(Rank.SPECIES), GSGoalKey.DB),
	/** RefSeq status values that restrict the considered genomic accessions. */
	@MDDescription("The refseq status values restrict the considered genomic accessions with respect to the given values. By default all values are allowed / included.")
	RES_SEQ_STATUS("refseq.status", new RefSeqStatusConfigInfo(Arrays.asList(RefSeqStatus.values())), GSGoalKey.DB),
	/** Whether extracted fasta files are gzipped. */
	@MDDescription("Whether to create gzipped extracted fasta files in goal `extractrefseqfasta`.")
	EXTRACT_REFSEQ_GZIP("reqseq.extract.gzip", new BooleanConfigParamInfo(false)),
	/** Whether fastq output of the filter and match goals is gzipped. */
	@MDDescription("Gzip fastq output from goals `filter` and `match`.")
	GZIP_FASTQ_OUTPUT("gzipFastqOutput", new BooleanConfigParamInfo(true), GSGoalKey.FILTER, GSGoalKey.MATCH, GSGoalKey.FASTA2FASTQ, GSGoalKey.DB2FASTQ),

	// Genbank data selection
	/** Maximum number of Genbank fasta files used per tax id. */
	@MDDescription("Determines the maximum number of fasta files used from Genbank per requested tax id. "
			+ "If this value is <= 0 then all fasta files will be used. "
			+ "Otherwise, if the corresponding number of matching files exceeds `genbank.maxPerTaxid`, then  best ones according to `genbank.fastaQualities` will be retained while adhering to this maximum.")
	MAX_FROM_GENBANK("genbank.maxPerTaxid", new IntConfigParamInfo(-1, Integer.MAX_VALUE, 1), GSGoalKey.DB),
	/** Allowed quality levels of Genbank fasta files. */
	@MDDescription("Determines the allowed quality levels of fasta files from Genbank. "
			+ "The values must be comma-separated. If a corresponding value is included in the list, "
			+ "then a fasta file for a requested tax id on that quality level will be included, "
			+ "otherwise not (while also respecting the conditions exerted via the keys `refSeq.limitForGenbankAccess` and `genbank.maxPerTaxid`). "
			+ "The quality levels are based on Genbank's [Assembly Summary File](https://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank.txt) (columns `version_status` and `assembly_level`). "
			+ "If the list is empty then no fasta files from Genbank will qualify.")
	FASTA_QUALITIES("genbank.fastaQualities", new FastaQualitiesConfigInfo(Arrays.asList(AssemblyQuality.COMPLETE_LATEST, AssemblyQuality.CHROMOSOME_LATEST)), GSGoalKey.DB),
	/** Whether only reference genomes are accepted. */
	@MDDescription("Whether only reference genomes are accepted or not. (Reference Genomes must be fetched from GenBank.)")
	REF_GEN_ONLY("genbank.referenceOnly", new BooleanConfigParamInfo(false), false, GSGoalKey.DB),

	// Database generation
	/** Maximum dust value above which low-complexity k-mers are dropped. */
	@MDDescription("When generating a database via the goal `db`, any low-complexity *k*-mer with too many "
			+ "repetitive sequences of base pairs may be omitted for storing. To do so, Genestrip employs a simple "
			+ "[genetic dust-filter](https://pubmed.ncbi.nlm.nih.gov/16796549/) for *k*-mers: "
			+ "It assigns a dust value *d* to each *k*-mer, and if *d* >  `maxDust`, then the *k*-mer will not be stored. "
			+ "Let *k(i)* be length of a *k*-mer's *i*-th substring s<sub>i</sub> of maximum length such that *s<sub>i</sub>(j) = s<sub>i</sub>(j-1)* holds for all bases in *s*. "
			+ "Given a *k*-mer with *n* such non-overlapping substrings and their lengths *k(1), ..., k(n)*, "
			+ "then *d = fib(k(1)) + ... + fib(k(n))*, where *fib(k(i))* is the Fibonacci number of *k(i)*. "
			+ "(The Fibonachi numbers are *fib(1) = 0*, *fib(2) = 1*, *fib(n) = fib(n-1) + fib(n-2)*.) "
			+ "E.g., for the *8*-mer `TTTCGCGA`, we have *n = 3* with *k(1) = 3* for `TTT`, *k(2) = 4* for `CGCG` and *k(3) = 1* for `A` which gives *d = fib(3) + fib(4) + fib(1) = 1 + 2 + 0 = 3*. "
			+ "For practical concerns `maxDust = 500` may be suitable. In this case, if *31*-mers were uniformly, randomly generated, "
			+ "then less than 0.00002 % of them would be dropped. If `maxDust = -1`, then dust-filtering is inactive.")
	MAX_DUST("maxDust", new IntConfigParamInfo(-1, Integer.MAX_VALUE, -1), GSGoalKey.DB),
	/** False positive probability of the temporary Bloom filter used by the tempindex goal. */
	TEMP_BLOOM_FILTER_FPP("tempBloomFilterFpp", new DoubleConfigParamInfo(0, 1, 0.001d, true), true, GSGoalKey.TEMPINDEX),
	/** Scaling factor applied to the estimated k-mer count when sizing the store. */
	@MDDescription("A scaling factor applied to the pre-computed *k*-mer count estimate (from the goal `fillsize`) to determine the allocated size of the *k*-mer store before filling it. "
			+ "A value greater than `1.0` reserves more space than the estimate; a value less than `1.0` reserves less. "
			+ "The default `1.0` uses the estimate as-is. Adjusting this value can be useful if the estimate from `fillsize` is slightly off.")
	// Exclusive lower bound: a factor of 0 would size the store to 0 -> a silently empty DB.
	DB_RESIZING_FACTOR("dbResizingFactor", new DoubleConfigParamInfo(0, Double.MAX_VALUE, 1, true), GSGoalKey.DB),
	/** False positive probability of the filtering database's Bloom filter. */
	@MDDescription("False positive probability (FPP) of the Bloom filter embedded in the filtering database (used by the goal `filter`). "
			+ "A lower value reduces false positives during filtering at the cost of a larger filter.")
	INDEX_BLOOM_FILTER_FPP("indexBloomFilterFpp", new DoubleConfigParamInfo(0, 1, GSConfigKey.INDEX_BLOOM_FILTER_FPP_DEFAULT, true), true, GSGoalKey.FILL_DB),
	/** False positive probability of the temporary Bloom filter used during the fill phase. */
	@MDDescription("False positive probability (FPP) of the temporary Bloom filter used inside the *k*-mer store during the database fill phase. "
			+ "This filter is discarded after sorting and replaced by the one governed by `optBloomFilterFpp`.")
	FILL_BLOOM_FILTER_FPP("fillBloomFilterFpp", new DoubleConfigParamInfo(0, 1, GSConfigKey.FILL_BLOOM_FILTER_FPP_DEFAULT, true), true, GSGoalKey.FILL_DB),
	/** False positive probability of the final matching database's Bloom filter. */
	@MDDescription("False positive probability (FPP) of the Bloom filter embedded in the final matching database after the *k*-mer store has been sorted and optimized. "
			+ "This is the filter used during matching when `useBloomFilterForMatch=true`.")
	OPT_BLOOM_FILTER_FPP("optBloomFilterFpp", new DoubleConfigParamInfo(0, 1, BlockedKMerBloomFilter.DEFAULT_FPP, true), true, GSGoalKey.FILL_DB),
	/** Whether the database uses the radix-indexed k-mer store. */
	@MDDescription("If `true`, the database's *k*-mer store uses the radix-indexed `RadixKMerStore` instead of the default sorted-array store. "
			+ "It is sized per radix bucket from the deduplicated per-bucket *k*-mer counts (see goal `tempindex`) and tends to be faster for lookups on large databases that exceed the CPU cache.")
	USE_RADIX_STORE("useRadixStore", new BooleanConfigParamInfo(true), false, GSGoalKey.FILL_DB),
	/** Number of low k-mer bits used as the radix index of the RadixKMerStore. */
	@MDDescription("Number of low *k*-mer bits used as the radix index of the `RadixKMerStore` (only relevant when `useRadixStore` is `true`). "
			+ "The store has `2^radixStoreBits` buckets; more bits give smaller, more cache-friendly buckets at the cost of a larger radix table. "
			+ "It also raises the store's value capacity (`MAX_VALUES`, the number of distinct values it can hold): each *k*-mer entry reserves `62 - radixStoreBits` low bits for the remaining *k*-mer bits (sized for the worst-case `k=31`) and uses the high `2 + radixStoreBits` bits (capped at 30) for the value index, so a wider radix leaves more bits for values. "
			+ "Because the value capacity grows with `radixStoreBits`, so does the memory of the store's value-index array; this scales with the larger databases that warrant a wider radix.")
	RADIX_STORE_BITS("radixStoreBits", new IntConfigParamInfo(RadixKMerStore.MIN_RADIX_BITS, 24, RadixKMerStore.DEFAULT_RADIX_BITS), false, GSGoalKey.TEMPINDEX, GSGoalKey.FILL_DB),
	/** Whether to XOR-combine the Bloom filter hash functions. */
	XOR_BLOOM_HASH("xorBloomHash", new BooleanConfigParamInfo(true)),
	/** Line length in bytes for generated fasta files. */
	FASTA_LINE_SIZE_BYTES("fastaLineSizeBytes", new IntConfigParamInfo(4096, 65536, 4096), true, GSGoalKey.DB),
	/** Whether the least-common-ancestor update uses only the database's own genomes. */
	@MDDescription("Perform database update regarding least common ancestors only based on genomes of tax ids as selected for the database generation (and not via all of a super-kingdom's RefSeq genomes).")
	MIN_UPDATE("minUpdate", new BooleanConfigParamInfo(false), false, GSGoalKey.UPDATE_DB),
	/** Whether to consider only complete-genome accession prefixes when updating the database. */
	@MDDescription("If `true`, then only genomic accessions with the prefixes `AC`, `NC_`, `NZ_` will be considered when updating the database. "
			+ "Otherwise, all genomic accessions will be considered for the update phase. See [RefSeq accession numbers and molecule types](https://www.ncbi.nlm.nih.gov/books/NBK21091/table/ch18.T.refseq_accession_numbers_and_mole/) for details.")
	UPDATE_WITH_COMPLETE_GENOMES_ONLY("refseq.updateWithCompleteGenomesOnly", new BooleanConfigParamInfo(false), GSGoalKey.UPDATE_DB),
	/** Whether to delete the temporary database after saving the final one. */
	@MDDescription("Whether to delete the temporary database after the final database has been saved or not.")
	REMOVE_TEMP_DB("removeTempDB", new BooleanConfigParamInfo(true), false, GSGoalKey.DB),
	/** Store only every stepSize-th k-mer of a genome. */
	@MDDescription("Stores *k*-mers in steps of `stepSize`. " +
			"E.g. if `stepSize=2` then only every second *k*-mer from a genome is considered for entry into the database.")
	STEP_SIZE("stepSize", new IntConfigParamInfo(1, Integer.MAX_VALUE, 1), GSGoalKey.DB),
	/** Whether to add artificial DATA-rank nodes as intermediate children of tax ids. */
	@MDDescription("Whether to add artificial nodes of rank `DATA` in the taxonomy tree as intermediate children of tax id nodes when filling the database. "
			+ "K-mers are then assigned to these intermediate DATA nodes rather than directly to the tax id, enabling finer-grained attribution within a taxon. "
			+ "Analogous to `idNodes` and `fileNodes`. BEWARE: This may cause a database build to fail as only up to 32767 tax ids are allowed.")
	DATA_NODES("dataNodes", new BooleanConfigParamInfo(false), false, GSGoalKey.DB),
	/** Whether to add artificial tax nodes for ids from fasta info lines. */
	@MDDescription("Whether to add artificial nodes in the tax tree to represent ids after '>' from fasta info lines for *k*-mers. BEWARE: This may cause a database build to fail as only up to 32767 tax ids are allowed.")
	ID_NODES("idNodes", new BooleanConfigParamInfo(false), false, GSGoalKey.DB),
	/** Whether to add artificial tax nodes representing fasta files. */
	@MDDescription("Whether to add artificial nodes in the tax tree to represent fasta files for *k*-mers.")
	FILE_NODES("fileNodes", new BooleanConfigParamInfo(false), false, GSGoalKey.DB),
	/** Whether lowercase bases are accepted for k-mers. */
	@MDDescription("Whether to accept lowercase bases for *k*-mers.")
	ENABLE_LOWERCASE_BASES("lowerCaseBases", new BooleanConfigParamInfo(true), false, GSGoalKey.DB),

	// SVG Database content
	/** Font name for texts in the generated tree. */
	@MDDescription("The font name for the texts in the generated tree.")
	SVG_FONT("svgFont", new StringConfigParamInfo("SansSerif"), false, GSGoalKey.SVG_TAX_TREE),
	/** Font size for texts in the generated tree. */
	@MDDescription("The font size for the texts in the generated tree.")
	SVG_FONT_SIZE("svgFontSize", new IntConfigParamInfo(1, 100, 18), false, GSGoalKey.SVG_TAX_TREE),
	/** Line-height factor for text in the generated tree. */
	@MDDescription("How much a line of text in the tree is de- or increased with regard to the normally required line height.")
	SVG_LINE_HEIGHT_FACTOR("svgLineHeightFactor", new DoubleConfigParamInfo(0.5, 10, 1), false, GSGoalKey.SVG_TAX_TREE),
	/** Indentation factor for child nodes in the tree. */
	@MDDescription("Factor for standard indentation of child nodes in tree.")
	SVG_INDENTATION_FACTOR("svgIndentFactor", new DoubleConfigParamInfo(0, 10, 0.75), false, GSGoalKey.SVG_TAX_TREE),
	/** Gap between a node's connector line and its text, relative to font size. */
	@MDDescription("Gap between horizontal line for child node and node text as a ratio of the font size.")
	SVG_TEXT_GAP_FACTOR("svgTextGapFactor", new DoubleConfigParamInfo(0, 1, 0.25), false, GSGoalKey.SVG_TAX_TREE),
	/** Additional indentation factor reflecting a node's k-mers. */
	@MDDescription("Factor for additional indentation to reflect *k*-mers of node. The base value is normalized to [0,1], where 1 corresponds to the maximum *k*-mers per taxid as stored in the database.")
	SVG_KMER_NODE_INDENT_FACTOR("svgKmerNodeIndentFactor", new DoubleConfigParamInfo(0, Double.MAX_VALUE, 0),false, GSGoalKey.SVG_TAX_TREE),
	/** Whether node indentation is based on evolutionary distance. */
	@MDDescription("Whether to perform indentation based on the evolutionary distance *d* for `svgKmerNodeIndentFactor` instead of the number of differing *k*-mers.")
	SVG_DISTANCE_INDENT("svgDistanceIndent", new BooleanConfigParamInfo(false), false, GSGoalKey.SVG_TAX_TREE),
	/** Whether requested tax ids are shown in bold. */
	@MDDescription("Whether to use bold text for tax ids requested via the project file `taxids.txt`.")
	SVG_REQ_NODES("svgReqNodesBold", new BooleanConfigParamInfo(true), false, GSGoalKey.SVG_TAX_TREE),
	/** Whether the rank is shown in the node text. */
	@MDDescription("Whether to add the rank in the node text.")
	SVG_SHOW_RANK("svgShowRank", new BooleanConfigParamInfo(false), false, GSGoalKey.SVG_TAX_TREE),
	/** Distance threshold above which a dashed line indicates unreliable measurement. */
	@MDDescription("Indicates a too large distance for reliable measuring via a dashed horizontal line")
	SVG_TOO_LARGE_DISTANCE("svgTooLargeDistance", new DoubleConfigParamInfo(0, 1, 1), false, GSGoalKey.SVG_TAX_TREE),
	/** Whether the longest-distance path branch is marked in red. */
	@MDDescription("Marks the branch of the longest distance path from a node, as used to compute distance-based indentation in red.")
	SVG_MARK_LONGEST_PATH("svgMarkLongestPath", new BooleanConfigParamInfo(false), false, GSGoalKey.SVG_TAX_TREE),
	/** Whether the distance value is printed after a node's name. */
	@MDDescription("Whether to print the distance value behind a node's name.")
	SVG_SHOW_DISTANCE("svgShowDistance", new BooleanConfigParamInfo(false), false, GSGoalKey.SVG_TAX_TREE),
	/** Whether the distance portion value is printed after a node's name. */
	@MDDescription("Whether to print the distance portion value behind a node's name.")
	SVG_SHOW_DISTANCE_PORTION("svgShowDistancePortion", new BooleanConfigParamInfo(false), false, GSGoalKey.SVG_TAX_TREE),


	// Match
	/** Number of reads per file between matching-progress log messages. */
	@MDDescription("Affects the log level `trace`: Defines after how many reads per fastq file, information on the matching progress is logged. If less than 1, then no progress information is logged.")
	LOG_PROGRESS_UPDATE_CYCLE("logProgressUpdateCycle", new LongConfigParamInfo(0, Long.MAX_VALUE, 1000000), GSGoalKey.MATCH, GSGoalKey.MATCHLR, GSGoalKey.FILTER),
	/** Whether to perform Kraken-style read classification. */
	@MDDescription("Whether to do read classification in the style of Kraken and KrakenUniq. Matching is faster without "
			+ "read classification and the columns `kmers`, `unique kmers` and `max contig length` in resulting CSV files are usually more conclusive anyways - "
			+ "in particular with respect to long reads. When read classification is off, the columns `reads` and `kmers from reads` will be 0 in resulting CSV files.")
	CLASSIFY_READS("classifyReads", new BooleanConfigParamInfo(true), GSGoalKey.MATCH),
	/** Whether unique k-mers are counted and reported. */
	@MDDescription("If `true`, the number of unique *k*-mers will be counted and reported. This requires less than 5% of additional main memory.")
	COUNT_UNIQUE_KMERS("countUniqueKMers", new BooleanConfigParamInfo(true), GSGoalKey.MATCH, GSGoalKey.MATCHLR),
	/** Whether the match goal also writes filtered fastq files. */
	@MDDescription("If `true`, then the goal `match` writes filtered fastq files in the same way that the goal `filter` does.")
	WRITE_FILTERED_FASTQ("writeFilteredFastq", new BooleanConfigParamInfo(false), GSGoalKey.MATCH, GSGoalKey.MATCHLR),
	/** Whether Kraken-style {@code .out} output files are written. */
	@MDDescription("If `true`, Genestrip will write output files with suffix `.out` in the [Kraken output format](https://ccb.jhu.edu/software/kraken/MANUAL.html#output-format) "
			+ "under `<base dir>/projects/<project_name>/krakenout` covering all reads with at least one matching *k*-mer.")
	WRITE_KRAKEN_STYLE_OUT("writeKrakenStyleOut", new BooleanConfigParamInfo(false), GSGoalKey.MATCH, GSGoalKey.MATCHLR),
	/** Whether all reads, not only classified ones, are written to Kraken-style output. */
	@MDDescription("If `false`, Genestrip will write only classified reads to kraken style output files.")
	WRITE_ALL("writeAll", new BooleanConfigParamInfo(true), GSGoalKey.MATCH),
	/** Whether a Bloom filter is used during matching. */
	@MDDescription("If `true` a bloom filter will be loaded and used during fastq file analysis (i.e. matching). "
			+ "Using the bloom filter tends to shorten matching time, if the most part of the reads cannot be classified because they contain *no* *k*-mers from the database. "
			+ "Otherwise, using the bloom filter might increase matching time by up to 30%. It also requires more main memory.")
	USE_BLOOM_FILTER_FOR_MATCH("useBloomFilterForMatch", new BooleanConfigParamInfo(true), GSGoalKey.MATCH, GSGoalKey.MATCHLR),
	/** Maximum count of a read's k-mers absent from the database for it to be classified. */
	@MDDescription("The absolute or relative maximum number of *k*-mers that do not need to be in the database for a read to be classified (read error count). "
			+ "If the number is above `maxReadTaxErrorCount`, then the read will not be classified. "
			+ "Otherwise the read will be classified in the same way as [done by Kraken](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2014-15-3-r46/figures/1). "
			+ "If `maxReadTaxErrorCount` is >= 1, then it is interpreted as an absolute number of *k*-mers. "
			+ "If >= 0 and < 1, it is interpreted as the ratio between the *k*-mers not in the database and all *k*-mers of the read. "
			+ "If `maxReadTaxErrorCount` < 0, then the read error count is disregarded, which means that even a single matching *k*-mer will lead to the read's classification.")
	MAX_READ_TAX_ERROR_COUNT("maxReadTaxErrorCount", new DoubleConfigParamInfo(-1, Double.MAX_VALUE, -1),
			GSGoalKey.MATCH, GSGoalKey.MATCHLR),
	/** Maximum count of a read's k-mers inconsistent with its class for it to be classified. */
	@MDDescription("The absolute or relative maximum number of *k*-mers that do not need to be consistent with a read's destined class for the read to be classified (read class error count). "
			+ "If the number is above `maxReadClassErrorCount`, then the read will not be classified. "
			+ "Otherwise the read will be classified in the same way as [done by Kraken](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2014-15-3-r46/figures/1). "
			+ "If `maxReadClassErrorCount` is >= 1, then it is interpreted as an absolute number of *k*-mers. "
			+ "If >= 0 and < 1, it is interpreted as the ratio between the inconsistent *k*-mers and all *k*-mers of the read. "
			+ "If `maxReadClassErrorCount` < 0, then the read error count is disregarded, which means that even a single matching *k*-mer will lead to the read's classification.")
	MAX_READ_CLASS_ERROR_COUNT("maxReadClassErrorCount", new DoubleConfigParamInfo(-1, Double.MAX_VALUE, -1),
			GSGoalKey.MATCH, GSGoalKey.MATCHLR),
	/** Minimum total of k-mers under a taxon required to classify a read to it. */
	@MDDescription("Can be set to adjust the minimal total of *k*-mers under taxon *t* required for a read to be classified to *t*. I.e., given a read *r* and taxon *t1* on the genus rank with two *k*-mers from *r* and taxon *t2* subordinate to *t1* on the species rank with one *k*-mer from *r*. Furthermore, *r* shall have no other *k*-mers matching any taxons. Then, if `minKMersForClass = 2`, *r* would not be classified to *t1* but to *t2* instead since the single *k*-mer under *t1* is below the threshold but the total of three *k*-mers under *t2* suffice.")
	MIN_KMERS_FOR_CLASS("minKMersForClass", new IntConfigParamInfo(1, Integer.MAX_VALUE, 1)),
	/** Size of the work queue between reader and consumer threads. */
	THREAD_QUEUE_SIZE("threadQueueSize", new IntConfigParamInfo(1, 10000, 1000), true),
	/** Initial read buffer size in bytes. */
	INITIAL_READ_SIZE_BYTES("initialReadSizeBytes", new IntConfigParamInfo(256, 65536, 4096), true),
	/** Maximum number of classification paths tracked per read. */
	MAX_CLASSIFICATION_PATHS("maxClassificationPaths", new IntConfigParamInfo(1, 128, 10), true),

	// Filter
	/** Whether the filter goal also writes a fastq file of dumped (unmatched) reads. */
	@MDDescription("If `true`, then `filter` will also generate a fastq file `dumped_...` with all reads not written to the corresponding filtered fastq file.")
	WRITE_DUMPED_FASTQ("writeDumpedFastq", new BooleanConfigParamInfo(false), GSGoalKey.FILTER),
	/** Minimum number of a read's k-mers in the Bloom index for it to pass the filter. */
	@MDDescription("The mininum number of a read's *k*-mers to be found in the bloom index such that the read is added to the filtered fastq file. If `minPosCountFilter=0`, then `posRatioFilter` becomes effective.")
	MIN_POS_COUNT_FILTER("minPosCountFilter", new IntConfigParamInfo(0, 1024, 1), GSGoalKey.FILTER),
	/** Minimum ratio of a read's k-mers in the Bloom index for it to pass the filter. */
	@MDDescription("Only effective if `minPosCountFilter=0`: The mininum ratio of a read's *k*-mers to be found in the bloom index such that the read is added to the filtered fastq file.")
	POS_RATIO_FILTER("posRatioFilter", new DoubleConfigParamInfo(0, 1, 0.2), GSGoalKey.FILTER),
	/** Whether base-pair probabilities are processed and written. */
	@MDDescription("Whether to process bp probabilities and potentially write them to filtered fastq files. (Takes slightly more memory if `true`.)")
	WITH_PROBS("withProbs", new BooleanConfigParamInfo(false),  GSGoalKey.FILTER, GSGoalKey.MATCH, GSGoalKey.MATCHLR),

	// DB 2 Fastq Goal
	/** List of tax ids, optionally suffixed with {@code +} to include descendants. */
	@MDDescription("List of tax ids separated by `,`. A tax id may have the suffix `+`, which means that taxonomic descendants from the project's database will be included. This list can alternatively be set via the command line parameter `-tx`.")
	TAX_IDS("taxids", new ListConfigParamInfo<String>(Collections.emptyList()) {
		@Override
		protected List<String> fromString(String s) {
			List<String> res = new ArrayList<>();
			if (s != null && !s.isEmpty()) {
				String[] ss = s.split(",");
				for (int i = 0; i < ss.length; i++) {
					String t = ss[i].trim();
					if (!t.isEmpty()) {
						res.add(t);
					}
				}
			}
			return res;
		}

		@Override
		public String getTypeDescriptor() {
			return "list of Strings";
		}
	}, GSGoalKey.DB2FASTQ_TAXIDS, GSGoalKey.DB2FASTQ),

	// Kraken
	/** Name of the Kraken/KrakenUniq binary. */
	KRAKEN_BIN("krakenBin", new StringConfigParamInfo("krakenuniq"), true),
	/** Name of the Kraken database. */
	KRAKEN_DB("krakenDB", new StringConfigParamInfo("krakenuniq"), true),
	/** Command-line expression template for invoking Kraken. */
	KRAKEN_EXEC_EXPR("krakenExecExpr", new StringConfigParamInfo("{0} -db {1} {2}"), true);

	/** Default false positive probability for the fill-phase Bloom filter. */
	public static final double FILL_BLOOM_FILTER_FPP_DEFAULT =  0.00000000001d;
	/** Default false positive probability for the filtering index Bloom filter. */
	public static final double INDEX_BLOOM_FILTER_FPP_DEFAULT = 0.00000001d;

	private final String name;
	private final ConfigParamInfo<?> param;
	private final boolean internal;
	private final GSGoalKey[] forGoals;

	GSConfigKey(String name, ConfigParamInfo<?> param, GSGoalKey... forGoals) {
		this(name, param, false, forGoals);
	}

	GSConfigKey(String name, ConfigParamInfo<?> param, boolean internal, GSGoalKey... forGoals) {
		if (name == null) {
			throw new NullPointerException("name");
		}
		if (param == null) {
			throw new NullPointerException("param");
		}
		this.name = name;
		this.param = param;
		this.internal = internal;
		this.forGoals = forGoals;
	}

	/**
	 * Indicates whether this configuration key is internal and hidden from generated documentation.
	 *
	 * @return whether this configuration key is internal and hidden from generated documentation
	 */
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

	/**
	 * @return whether this configuration key applies to the given goal; a {@code null} goal or a key
	 *         declared for no specific goal matches any goal
	 */
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

	/**
	 * Prints a Markdown table describing all non-internal configuration keys that apply to the given
	 * goal (or all such keys if {@code filterGoalKey} is {@code null}).
	 *
	 * @param ps the stream to print the Markdown table to
	 * @param filterGoalKey the goal to filter keys by, or {@code null} for all keys
	 */
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

	/**
	 * The kind of RefSeq sequence files to include (genomic, RNA variants or all).
	 */
	public enum SeqType {
		/** Genomic sequences. */
		GENOMIC,
		/** RNA sequences. */
		RNA,
		/** Messenger RNA (mRNA) sequences. */
		M_RNA,
		/** All RNA sequence variants. */
		ALL_RNA,
		/** All sequence types. */
		ALL;

		/**
		 * Returns the sequence type with the given name.
		 *
		 * @param name the name to look up (case-insensitive)
		 * @return the sequence type whose name matches (case-insensitively), or {@code null} if none
		 */
		public static SeqType byName(String name) {
			for (SeqType r : SeqType.values()) {
				if (r.name().equalsIgnoreCase(name)) {
					return r;
				}
			}
			return null;
		}
	}

	/**
	 * RefSeq status values used to restrict which genomic accessions are considered.
	 */
	public enum RefSeqStatus {
		/** Status not available. */
		NA("na"),
		/** Unknown status. */
		UNKNOWN("UNKNOWN"),
		/** Reviewed status. */
		REVIEWED("REVIEWED"),
		/** Validated status. */
		VALIDATED("REVIEWED"),
		/** Provisional status. */
		PROVISIONAL("REVIEWED"),
		/** Predicted status. */
		PREDICTED("REVIEWED"),
		/** Inferred status. */
		INFERRED("REVIEWED"),
		/** Model status. */
		MODEL("REVIEWED");

		private final String name;

		private RefSeqStatus(String name) {
			this.name = name;
		}

		/**
		 * Returns the RefSeq status name associated with this status value.
		 *
		 * @return the RefSeq status name
		 */
		public String getName() {
			return name;
		}
	}

	/**
	 * Parameter info for a comma/semicolon-separated list of Genbank {@link AssemblyQuality} values.
	 */
	public static class FastaQualitiesConfigInfo extends ListConfigParamInfo<AssemblyQuality> {
		/**
		 * Creates parameter info for a list of Genbank quality values.
		 *
		 * @param defaults the default list of quality values
		 */
		public FastaQualitiesConfigInfo(List<AssemblyQuality> defaults) {
			super(defaults);
		}

		@Override
		protected List<AssemblyQuality> fromString(String qs) {
			List<AssemblyQuality> res = new ArrayList<>();
			if (qs != null) {
				StringTokenizer tokenizer = new StringTokenizer(qs, ",;");
				while (tokenizer.hasMoreTokens()) {
					try {
						res.add(AssemblyQuality.valueOf(tokenizer.nextToken().trim()));
					} catch (IllegalArgumentException e) {
						// Unknown token: report as an invalid value (the framework falls back to default).
						return null;
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

	/**
	 * Parameter info for a log level, restricting the value to the recognized level names.
	 */
	public static class LogLevelConfigParamInfo extends StringConfigParamInfo {
		private final String[] LEVELS = new String[] { "all", "trace", "debug", "info", "warn", "error", "fatal", "off" };

		/**
		 * Creates parameter info for a log level.
		 *
		 * @param defaultValue the default log level
		 */
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

	/**
	 * Parameter info for a taxonomic {@link Rank} value (or {@code null} for no rank).
	 */
	public static class RankConfigParamInfo extends ConfigParamInfo<Rank> {
		/**
		 * Creates parameter info for a taxonomic rank.
		 *
		 * @param defaultValue the default rank, or {@code null} for none
		 */
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

	/**
	 * Parameter info for a {@link SeqType} value.
	 */
	public static class SeqTypeConfigParamInfo extends ConfigParamInfo<SeqType> {
		/**
		 * Creates parameter info for a sequence type.
		 *
		 * @param defaultValue the default sequence type
		 */
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

	/**
	 * Parameter info for a comma/semicolon-separated list of {@link RefSeqStatus} values.
	 */
	public static class RefSeqStatusConfigInfo extends ListConfigParamInfo<RefSeqStatus> {
		/**
		 * Creates parameter info for a list of RefSeq status values.
		 *
		 * @param defaults the default list of status values
		 */
		public RefSeqStatusConfigInfo(List<RefSeqStatus> defaults) {
			super(defaults);
		}

		@Override
		protected List<RefSeqStatus> fromString(String qs) {
			List<RefSeqStatus> res = new ArrayList<>();
			if (qs != null) {
				StringTokenizer tokenizer = new StringTokenizer(qs, ",;");
				while (tokenizer.hasMoreTokens()) {
					try {
						res.add(RefSeqStatus.valueOf(tokenizer.nextToken().trim()));
					} catch (IllegalArgumentException e) {
						// Unknown token: report as an invalid value (the framework falls back to default).
						return null;
					}
				}
			}
			return res;
		}

		@Override
		public String getMDRangeDescriptor() {
			StringBuilder builder = new StringBuilder();
			RefSeqStatus[] types = RefSeqStatus.values();
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
}
