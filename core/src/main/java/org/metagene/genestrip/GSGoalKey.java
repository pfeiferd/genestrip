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

import org.metagene.genestrip.make.MDDescription;
import org.metagene.genestrip.make.GoalKey;

/**
 * Enumeration of all Genestrip goals. Each constant binds a goal name and flags for whether the goal
 * is user-facing and whether it is included in transitive cleaning. Constants annotated with
 * {@link MDDescription} contribute to the generated documentation.
 */
public enum GSGoalKey implements GoalKey {
	/** Shows the user-related goals. */
	@MDDescription("Show user-related goals.")
	SHOW("show", true),
	/** Shows the user-related and most internal goals. */
	@MDDescription("Show user-related and most internal goals.")
	SHOWALL("showall", true),
	/** Generates the k-mer matching database and the filtering database for the given project. */
	@MDDescription("Generate the *k*-mer matching database and also the filtering database with respect to the given project.")
	GENALL("genall", true),
	/** Clears a project's {@code csv}, {@code log} and {@code krakenout} folders, deleting all their files. */
	@MDDescription("Clear the folders `csv`, `log` and `krakenout`  of a project. This will delete all files the respective folders!")
	CLEAR("clear", true),
	/** Generates the k-mer matching database for the given project. */
	@MDDescription("Generate the database for *k*-mer matching with respect to the given project.")
	DB("db", true),
	/** Writes information on a project's database content to a CSV file. */
	@MDDescription("Write information on a project's database content to a CSV file.")
	DBINFO("dbinfo", true),
	/** Loads config information for a database without necessarily loading the database itself. */
	@MDDescription("Load config information for a database - but not (necessarily) that database itself.")
	DBCONF("dbconf", false),
	/** Shows config information for a database without loading the database itself. */
	@MDDescription("Show config information for a database - without (necessarily) loading the database itself.")
	SHOWDBCONF("showdbconf", false),
	/** Writes information about a project's temporary database content to a CSV file. */
	@MDDescription("Write information about a project's temporary database content to a CSV file.")
	TEMP_DBINFO("tempdbinfo", false),
	/** Generates fastq files from the database, one per tax id, with each k-mer represented by a read. */
	@MDDescription("Generate fastq files from the database. A respective fastq file will contain all *k*-mers specifically associated with a "
			+ "single tax id from the database where each *k*-mer is represented by a read consisting of *k* bases. Respective fastq files will be stored "
			+ "in `<base dir>/projects/<project_name>/fastq` with the file name format `<project_name>_db2fastq_<taxid>.fastq.gz`. "
			+ "The command line option `tx` serves at selecting the corresponding tax ids for the fastq files to be generated. "
			+ "If the option is omitted, then fastq files for *all* tax ids from the database will be generated.")
	DB2FASTQ("db2fastq", true),
	/** Generates a filtering database for the given project. */
	@MDDescription("Generate a filtering database with respect to a given project.")
	INDEX("index", true),
	/** Analyzes fastq files given via the {@code -f} or {@code -m} option. */
	@MDDescription("Analyze fastq files as given by the `-f` or `-m` option. The resulting CSV file(s) will be stored in `<base dir>/projects/<project_name>/csv` unless specified otherwise via the `-r` option.")
	MATCH("match", true),
	/** Same as {@code match} but without read classification. */
	@MDDescription("Same as `match` but without doing read classification. This corresponds to the configuration setting `classifyReads=false`.")
	MATCHLR("matchlr", true),
	/** Filters fastq files given via the {@code -f} or {@code -m} option. */
	@MDDescription("Filter fastq files as given by the `-f` or `-m` option. The resulting filtered fastq file(s) `filtered_...` will be stored under `<base dir>/projects/<project_name>/fastq/` unless specified otherwise via the `-r` option.")
	FILTER("filter", true),
	/** Extracts reads from fastq files based on matching descriptors. */
	@MDDescription("Extract reads from fastq files based on matching descriptors. See also config key `extractKey`.")
	EXTRACT("extract", true),
	/** Generates an SVG visualization of the database's taxonomy tree. */
	@MDDescription("Generates a compact visual taxonomy tree in SVG-format to represent the database content.")
	SVG_TAX_TREE("svgtaxtree", true),

	// Non user Goals
	/** Transforms a fasta file or streaming resource to a fastq file. */
	@MDDescription("Transform a fasta file or streaming resource to fastq file.")
	FASTA2FASTQ("fasta2fastq"),
	/** Analyzes fastq files given via the {@code -f} or {@code -m} option. */
	@MDDescription("Analyze fastq files as given by the `-f` or `-m` option.")
	MATCHRES("matchres"),
	/** Same as {@code matchres} but without read classification. */
	@MDDescription("Same as `matchres` but without doing read classification.")
	MATCHRESLR("matchreslr"),
	/** Creates the data folders under {@code <base dir>/common}, including {@code common} itself. */
	@MDDescription("Create data folders in `<base dir>/common` including `common` itself.")
	COMMON_SETUP("commonsetup", false, false),
	/** Creates the data folders under {@code <base dir>/<project>}. */
	@MDDescription("Create data folders in `<base dir>/<project>`.")
	SETUP("setup", false, false),
	/** Downloads the taxonomy. */
	@MDDescription("Download the taxonomy.")
	TAXDOWNLOAD("taxdownload", false, false),
	/** Loads the taxonomy into memory. */
	@MDDescription("Load the taxonomy into memory.")
	TAXTREE("taxtree"),
	/** Computes the tax ids for the project's database. */
	@MDDescription("Compute the taxids for the project's database.")
	TAXNODES("taxnodes"),
	/** Downloads the RefSeq release number. */
	@MDDescription("Download the RefSeq release number.")
	REFSEQRELEASE("refseqrelease"),
	/** Stores the RefSeq release number in a properties file. */
	@MDDescription("Store the RefSeq release number in a properties file.")
	REFSEQ_PROP("refseqprop"),
	/** Downloads the RefSeq release catalog files. */
	@MDDescription("Download the RefSeq release catalog files.")
	REFSEQCAT("refseqcat", false, false),
	/** Loads MD5 checksums for RefSeq files into memory. */
	@MDDescription("Load MD5 checksums for RefSeq files into memory")
	CHECKSUMMAP("checksummap", false, false),
	/** Loads the RefSeq category names requested for the project's database. */
	@MDDescription("Load the RefSeq category names as requested for the project's database.")
	CATEGORIES("categories"),
	/** Downloads the genomic files from RefSeq for the requested categories. */
	@MDDescription("Download the genomic files from RefSeq for the requested categories.")
	REFSEQFNA("refseqfna", false, false),
	/** Computes the number of RefSeq accession entries to keep in memory. */
	@MDDescription("Compute the number of RefSeq accession entries to be kept in memory.")
	ACCMAPSIZE("accmapsize"),
	/** Loads the required RefSeq accession entries into memory. */
	@MDDescription("Load the required RefSeq accession entries into memory.")
	ACCMAP("accmap"),
	/** Determines the tax ids for which additional fasta files should be downloaded from Genbank. */
	@MDDescription("Determine for tax ids for which additional fasta files from Genbank should be downloaded.")
	TAXFROMGENBANK("taxfromgenbank"),
	/** Downloads Genbank's assembly catalog file. */
	@MDDescription("Download Genbank's assembly catalog file.")
	ASSEMBLYDOWNLOAD("assemblydownload", false, false),
	/** Determines the fasta files to be downloaded from Genbank. */
	@MDDescription("Determine the fasta files to be downloaded from Genbank.")
	FASTAGSENBANK("fastasgenbank"),
	/** Downloads the requested fasta files from Genbank. */
	@MDDescription("Download the requested fasta files from Genbank.")
	FASTAGSENBANKDL("fastasgenbankdl", false, false),
	/** Downloads additionally requested fasta files. */
	@MDDescription("Download additionally requested fasta files.")
	ADD_DOWNLOADS("adddownloads", false, false),
	/** Associates additional fasta files with tax ids in memory. */
	@MDDescription("Associate additional fasta files with tax ids in memory.")
	ADD_FASTAS("addfastas"),
	/** Precomputes the number of k-mers for the project's database. */
	@MDDescription("Precompute the number of *k*-mers for the project's database.")
	FILLSIZE("fillsize"),
	/** Fills the temporary bloom index with k-mers. */
	@MDDescription("Fill the temporary bloom index with *k*-mers.")
	TEMPINDEX("tempindex"),
	/** Fills the database with k-mers. */
	@MDDescription("Fill the database with *k*-mers.")
	FILL_DB("filldb"),
	/** Saves the temporary database. */
	@MDDescription("Save temporary database.")
	TEMPDB("tempdb"),
	/** Loads the temporary database into memory. */
	@MDDescription("Load the temporary database into memory.")
	LOAD_TEMPDB("loadtempdb"),
	/** Updates the temporary database with least common ancestors of k-mers. */
	@MDDescription("Update the temporary database with regard to least common ancestors of *k*-mers.")
	UPDATE_DB("updatedb"),
	/** Loads the matching database into memory. */
	@MDDescription("Load the matching database into memory.")
	LOAD_DB("loaddb"),
	/** Fills the filtering database with k-mers. */
	@MDDescription("Fill the filtering database with *k*-mers.")
	FILL_INDEX("fillindex"),
	/** Loads the filtering database into memory. */
	@MDDescription("Load the filtering database into memory.")
	LOAD_INDEX("loadindex"),
	/** Determines the tax ids for the {@code db2fastq} goal from the command line argument. */
	@MDDescription("Determine tax ids for the goal `db2fastq` from the command line argument.")
	DB2FASTQ_TAXIDS("db2fastqtaxids"),
	/** Determines the fastq mapping from command line arguments. */
	@MDDescription("Determine the fastq mapping from command line arguments.")
	FASTQ_MAP("fastqmap"),
	/** Transforms URLs of fastq files to be downloaded to local paths. */
	@MDDescription("Transform URLs of fastq files to be downloaded to local paths.")
	FASTQ_MAP_TRANSFORM("fastqmaptransform"),
	/** Determines the fasta mapping from command line arguments. */
	@MDDescription("Determine the fasta mapping from command line arguments.")
	FASTA_MAP("fastamap"),
	/** Transforms URLs of fasta files to be downloaded to local paths. */
	@MDDescription("Transform URLs of fasta files to be downloaded to local paths.")
	FASTA_MAP_TRANSFORM("fastamaptransform"),
	/** Downloads requested fastq files given via URLs. */
	@MDDescription("Download fastq files given via URLs as requested.")
	FASTQ_DOWNLOAD("fastqdownload", false, false),
	/** Downloads requested fasta files given via URLs. */
	@MDDescription("Download fasta files given via URLs as requested.")
	FASTA_DOWNLOAD("fastadownload", false, false),
	/** For internal use: invokes kraken and counts results. */
	@MDDescription("For internal use (to invoke kraken and count results).")
	KRAKENCOUNT("krakencount"),
	/** For internal use: writes kraken results to a file. */
	@MDDescription("For internal use (to write kraken results to a file).")
	KRAKENRES("krakenres"),
	/** Downloads and installs a project's database from a given URL. */
	@MDDescription("Download and install a project's database via a given URL.")
	DB_DOWNLOAD("dbdownload", false, false),
	/** Checks whether the downloaded RefSeq release matches the current release on the download server. */
	@MDDescription("Check whether the downloaded RefSeq release is equal to the current release on the download server.")
	CHECK_REFSEQ_RNUM("checkrefseqrnum"),
	/** Extracts gzipped fasta files individually from RefSeq's bundled offering for database entry. */
	@MDDescription("Extract gzipped fasta files as individual files from RefSeq's bundled offering as used for database entry.")
	EXTRACT_REFSEQ_FASTA("extractrefseqfasta"),
	/** Runs {@code extractrefseqfasta} and generates the matching CSV file. */
	@MDDescription("Run `extractrefseqfasta` and generate matching CSV file.")
	EXTRACT_REFSEQ_CSV("extractrefseqcsv");

	private final boolean forUser;
	private final String name;
	private final boolean transClean;

	private GSGoalKey(String name) {
		this(name, false, true);
	}

	private GSGoalKey(String name, boolean forUser) {
		this(name, forUser, true);
	}

	private GSGoalKey(String name, boolean forUser, boolean transClean) {
		this.name = name;
		this.forUser = forUser;
		this.transClean = transClean;
	}

	@Override
	public boolean isTransClean() {
		return transClean;
	}

	/**
	 * Returns whether this goal is intended to be invoked directly by users.
	 *
	 * @return whether this goal is intended to be invoked directly by users
	 */
	public boolean isForUser() {
		return forUser;
	}

	@Override
	public String getName() {
		return name;
	}

	@Override
	public String toString() {
		return name;
	}

	/**
	 * Prints a Markdown table listing all goals with their user/clean flags and descriptions.
	 *
	 * @param ps the print stream to write the Markdown table to
	 */
	public static void printGoalInfo(PrintStream ps) {
		ps.print('|');
		ps.print("Name");
		ps.print('|');
		ps.print("User Goal");
		ps.print('|');
		ps.print("Rec. Clean");
		ps.print('|');
		ps.print("Description");
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
		ps.println();

		for (GSGoalKey goalKey : GSGoalKey.values()) {
			ps.print('|');
			ps.print('`');
			ps.print(goalKey.getName());
			ps.print('`');
			ps.print('|');
			ps.print(goalKey.isForUser() ? "X" : "");
			ps.print('|');
			ps.print(goalKey.isTransClean() ? "X" : "");
			ps.print('|');
			Annotation[] annotations;
			try {
				annotations = GSGoalKey.class.getField(goalKey.name()).getAnnotations();
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
			ps.println();
		}
	}
}