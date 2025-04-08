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

import org.metagene.genestrip.goals.MDDescription;
import org.metagene.genestrip.make.GoalKey;

public enum GSGoalKey implements GoalKey {
	@MDDescription("Show user-related goals.")
	SHOW("show", true), 
	@MDDescription("Show user-related and most internal goals.")
	SHOWALL("showall", true),
	@MDDescription("Generate the *k*-mer matching database and also the filtering database with respect to the given project.")
	GENALL("genall", true),
	@MDDescription("Clear the folders `csv`, `db` and `krakenout`  of a project. This will delete all files the respective folders!")
	CLEAR("clear", true), 
	@MDDescription("Generate the database for *k*-mer matching with respect to the given project.")
	DB("db", true), 
	@MDDescription("Write information about a project's database content to a CSV file.")
	DBINFO("dbinfo", true),
	@MDDescription("Write information about a project's temporary database content to a CSV file.")
	TEMP_DBINFO("tempdbinfo", false),
	@MDDescription("Generate fastq files from the database. A respective fastq file will contain all *k*-mers specifically associated with a "
			+ "single tax id from the database where each *k*-mer is represented by a read consisting of *k* bases. Respective fastq files will be stored "
			+ "in `<base dir>/projects/<project_name>/fastq` with the file name format `<project_name>_db2fastq_<taxid>.fastq.gz`. "
			+ "The command line option `tx` serves at selecting the corresponding tax ids for the fastq files to be generated. "
			+ "If the option is omitted, then fastq files for *all* tax ids from the database will be generated.")
	DB2FASTQ("db2fastq", true),
	@MDDescription("Generate a filtering database with respect to a given project.")
	INDEX("index", true),
	@MDDescription("Analyze fastq files as given by the `-f` or `-m` option. The resulting CSV file(s) will be stored in `<base dir>/projects/<project_name>/csv` unless specified otherwise via the `-r` option.")
	MATCH("match", true),
	@MDDescription("Same as `match` but without doing read classification. This corresponds to the configuration setting `classifyReads=false`.")
	MATCHLR("matchlr", true),
	@MDDescription("Filter fastq files as given by the `-f` or `-m` option. The resulting filtered fastq file(s) `filtered_...` will be stored under `<base dir>/projects/<project_name>/fastq/` unless specified otherwise via the `-r` option.")
	FILTER("filter", true),

	// Non user Goals
	@MDDescription("Transform a fasta file or streaming resource to fastq file.")
	FASTA2FASTQ("fasta2fastq"),
	@MDDescription("Analyze fastq files as given by the `-f` or `-m` option.")
	MATCHRES("matchres"),
	@MDDescription("Same as `matchres` but without doing read classification.")
	MATCHRESLR("matchreslr"),
	@MDDescription("Create data folders in `<base dir>/common` including `common` itself.")
	COMMON_SETUP("commonsetup"), @MDDescription("Create data folders in `<base dir>/<project>`.")
	SETUP("setup"), @MDDescription("Download the taxonomy.")
	TAXDOWNLOAD("taxdownload"), @MDDescription("Load the taxonomy into memory.")
	TAXTREE("taxtree"), @MDDescription("Compute the taxids for the project's database.")
	TAXNODES("taxnodes"), @MDDescription("Download the RefSeq release number.")
	REFSEQRELEASE("refseqrelease"), @MDDescription("Download the RefSeq release catalog files.")
	REFSEQCAT("refseqcat"), @MDDescription("Load MD5 checksums for RefSeq files into memory")
	CHECKSUMMAP("checksummap"),
	@MDDescription("Load the RefSeq category names as requested for the project's database.")
	CATEGORIES("categories"), @MDDescription("Download the genomic files from RefSeq for the requested categories.")
	REFSEQFNA("refseqfna"), @MDDescription("Compute the number of RefSeq accession entries to be kept in memory.")
	ACCMAPSIZE("accmapsize"), @MDDescription("Load the required RefSeq accession entries into memory.")
	ACCMAP("accmap"),
	@MDDescription("Determine for tax ids for which additional fasta files from Genbank should be downloaded.")
	TAXFROMGENBANK("taxfromgenbank"), @MDDescription("Download Genbank's assembly catalog file.")
	ASSEMBLYDOWNLOAD("assemblydownload"), @MDDescription("Determine the fasta files to be downloaded from Genbank.")
	FASTAGSENBANK("fastasgenbank"), @MDDescription("Download the requested fasta files from Genbank.")
	FASTAGSENBANKDL("fastasgenbankdl"), @MDDescription("Download additionally requested fasta files.")
	ADD_DOWNLOADS("adddownloads"), @MDDescription("Associate additional fasta files with tax ids in memory.")
	ADD_FASTAS("addfastas"), @MDDescription("Precompute the number of *k*-mers for the project's database.")
	FILLSIZE("fillsize"), 
	@MDDescription("Fill the temporary bloom index with *k*-mers.")
	TEMPINDEX("tempindex"), 
	@MDDescription("Fill the database with *k*-mers.")
	FILL_DB("filldb"),
	@MDDescription("Save temporary database.")
	TEMPDB("tempdb"), 
	@MDDescription("Load the temporary database into memory.")
	LOAD_TEMPDB("loadtempdb"), 
	@MDDescription("Update the temporary database with regard to least commonn ancestors of *k*-mers.")
	UPDATE_DB("updatedb"), 
	@MDDescription("Load the matching database into memory.")
	LOAD_DB("loaddb"),
	@MDDescription("Fill the filtering database with *k*-mers.")
	FILL_INDEX("fillindex"), 
	@MDDescription("Load the filtering database into memory.")
	LOAD_INDEX("loadindex"),
	@MDDescription("Determine tax ids for the goal `db2fast` from the command line argument.")
	DB2FASTQ_TAXIDS("db2fastqtaxids"),
	@MDDescription("Determine the fastq mapping from command line arguments.")
	FASTQ_MAP("fastqmap"),
	@MDDescription("Transform URLs of fastq files to be downloaded to local paths.")
	FASTQ_MAP_TRANSFORM("fastqmaptransform"),
	@MDDescription("Determine the fasta mapping from command line arguments.")
	FASTA_MAP("fastamap"),
	@MDDescription("Transform URLs of fasta files to be downloaded to local paths.")
	FASTA_MAP_TRANSFORM("fastamaptransform"),
	@MDDescription("Download fastq files given via URLs as requested.")
	FASTQ_DOWNLOAD("fastqdownload"),
	@MDDescription("Download fasta files given via URLs as requested.")
	FASTA_DOWNLOAD("fastadownload"),
	@MDDescription("For internal use (to invoke kraken).")
	KRAKENRES("krakenres"),
	@MDDescription("Download and install a project's database via a given URL.")
	DB_DOWNLOAD("dbdownload");

	private final boolean forUser;
	private final String name;

	private GSGoalKey(String name) {
		this(name, false);
	}

	private GSGoalKey(String name, boolean forUser) {
		this.name = name;
		this.forUser = forUser;
	}

	public boolean isForUser() {
		return forUser;
	}

	public String getName() {
		return name;
	}

	public static void printGoalInfo(PrintStream ps) {
		ps.print('|');
		ps.print("Name");
		ps.print('|');
		ps.print("User Goal");
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
		ps.println();

		for (GSGoalKey goalKey : GSGoalKey.values()) {
			ps.print('|');
			ps.print('`');
			ps.print(goalKey.getName());
			ps.print('`');
			ps.print('|');
			ps.print(goalKey.isForUser() ? "X" : "");
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