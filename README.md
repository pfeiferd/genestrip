Genestrip - Hyper-efficient unique *k*-mer counting and sequence filtering for selected groups of species
===============================================

# Introduction

Metagenomic analysis has become an extremely import field in bio informatics. To analyse large sets of reads, researchers use highly efficient software tools based on *k*-mer matching and counting such as [Kraken](https://github.com/DerrickWood/kraken), [Kraken 2](https://github.com/DerrickWood/kraken2), [KrakenUniq](https://github.com/fbreitwieser/krakenuniq) or [Kaiju](https://github.com/bioinformatics-centre/kaiju). With regard to pathogen detection [KrakenUniq](https://github.com/fbreitwieser/krakenuniq) is particularly suited because of its
 [very low false positive rate](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-018-1568-0) but also because its efficiency, precision and sensitivity.

To avoid false positive classifications of reads or mismatches of *k*-mers, all of these tools resort to very large databases, containing millions of encoded *k*-mers along with their respective tax ids. A database is usually loaded entirely into main memory and consumes tens of gigabytes of space. This requirement mandates specialized and expensive compute servers with very large main memory.

If only a small group of a few dozen species or strains is to be considered, one could simply *generate a smaller database that contains just the k-mers of the respective genomes*. **However, when done improperly, this may lead to many false positives at matching time.** The reason for this is that due to the evolutionary relationship of organisms, *a large fraction of k-mers that belong to the genome of a species are unspecific* in a sense that they can as well be found in the genome other species not considered for analysis.

Genestrip offers an efficient and sophisticated database generation process that accounts for this problem: 
* Genstrips's databases can be generated in a fully automated way on the basis of the [RefSeq](https://ftp.ncbi.nlm.nih.gov/refseq/release/). 
* Species or strains, whose *k*-mers should be contained, are specified via a text file containing the corresponding [tax ids](https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/). 
* After the database generation has completed, it can be used for highly efficient *k*-mer matching.
* As a Genestrip database comprises just the genomes of the specified tax ids, it tends to remain small. Therefore analysis can be performed on edge machines as well as on small sized PCs.
* Much like KrakenUniq, Genestrip supports (exact) unique *k*-mer counting and total *k*-mer counting. Related results are stored in a CSV file. 
* Moreover, a fastq file can be *filtered* via the *k*-mers of a database to reduce the size of the file considerably. Still, all reads containing *k*-mers of the selected tax ids will be kept in the filtered fastq file. A similar functionality can be achieved via [biobloom](https://github.com/bcgsc/biobloom), but preparing and applying related filters is very convenient with Genestrip.

# License

[Genestrip is free for non-commercial use.](./LICENSE.txt) Please contact [daniel.pfeifer@progotec.de](mailto:daniel.pfeifer@progotec.de) if you are interested in a commercial license.

# Building and installing

Genestrip is structured as a standard [Maven 2 or 3](https://maven.apache.org/) project and is compatible with the [JDK](https://jdk.java.net/) 1.8 or higher.[^1]

To build it, `cd` to the installation directory `genestrip`. Given a matching Maven and JDK installation, `mvn install` will compile, test and install the Genestrip program library. Afterwards a self-contained and executable `genestrip.jar` file will be stored under `./lib`. 

[^1]: Counter to common belief, Java can well be used for such high performance applications when using its programming facilities the right way.

# Generating the sample database

The Genestrip installation holds additional folders that include the sample project `human_virus`. After building Genestrip, you may call
`sh ./bin/genestrip.sh human_virus dbinfo`
in order to generate the `human_virus` database and create a CSV file with basic information about the database content.

Genestrip follows a goal-oriented approach in order to create any result files (in similarity to [make](https://pubs.opengroup.org/onlinepubs/9699919799/utilities/make.html)). So, when generating the `human_virus` database, Genestrip will
1. download the [taxonomy](https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip) and unzip it to `./data/common`,
1. download the [refseq release catalog](https://ftp.ncbi.nlm.nih.gov/refseq/release/release-catalog/) to `./data/common/refseq`,
1. download [all virus related RefSeq fna files](https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/) to `./data/commmon/refseq` (which is currently just a single file),
1. perform several follow-up goals, until the database file `human_virus_db.ser.gz` is finally made under `./data/projects/human_virus/db`, 
1. create a CSV file `human_virus_dbinfo.csv` under `./data/projects/human_virus/csv`, which contains basic information about the database, i.e. the number of *k*-mers stored per tax id.

The generated database comprises *k*-mers for all viruses according to the tax id file `./data/project/human_virus/taxids.txt`.


# Generating your own database

Generating your own database is straight-forward:
1. Create a folder `<project_name>` under `./data/projects/`. This is the place for all subsequent files that specify the content of a database to be generated. It is also the core name of a related database.
1. Create a text file `./data/projects/<project_name>/taxids.txt` with one tax id per line. The database will *only* contain *k*-mers from genomes that belong to tax ids referenced in the file `taxids.txt`. If `taxids.txt` contains a tax id from a non-leaf node of the [taxonomy](https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip), then all subordinate tax ids and *k*-mers from respective [RefSeq](https://ftp.ncbi.nlm.nih.gov/refseq/release/) genomes will be included in the database as well.
1. Create a text file `./data/projects/<project_name>/categories.txt`. This file tells Genestrip, which categories of organisms should be *considered* when generating the database and when determining a least common ancestor tax id for *k*-mer in the finalized database. You also have to ensure that the categories from `categories.txt` comprise all of the tax ids from `taxids.txt`: Only genomes of comprised tax ids will be used to generate the database.
1. The following categories are allowed, and originate from the [RefSeq](https://ftp.ncbi.nlm.nih.gov/refseq/release/):
`archaea`, `bacteria`, `complete`, `fungi`, `invertebrate`, `mitochondrion`, `other`, `plant`, `plasmid`, `plastid`, `protozoa`, `vertebrate_mammalian`, `vertebrate_other` and `viral`. You may enter one category per line.
The *entire* set of genomes that belong a category referenced in `categories.txt` will be used to determine the least common ancestors tax ids of the *k*-mers stored in the finalized database. The more categories you reference in `categories.txt`, the more genomes files will have to be downloaded from the [RefSeq](https://ftp.ncbi.nlm.nih.gov/refseq/release/) and the longer the database generation process will take. E.g., if only the category `viral` is included, it should just take a few minutes, but with `bacteria` included, it may typically take about a day or more.
So you should chose a suitable set of categories to balance the completeness of considered genomes and efficiency of the database generation process.

1. Start the database generation process via the command `sh ./bin/genestrip.sh <project_name> dbinfo`. This will also create a CSV file `./data/projects/<project_name>/csv/<project_name>_dbinfo.csv` with basic information about the database, i.e. the number of *k*-mers stored per tax id. The path of the database will be `./data/projects/<project_name>/db/<project_name>_db.ser.gz`

In general, Genestrip organizes a project folder `./data/projects/<project_name>` by means of the following sub-folders:
* `csv` is where analysis results of fastq files will be stored (by default).
* `fastq` is where Genestrip will store filtered fastq files. You may also put the fastq files to be analyzed there (but they can be read from any other path too).
* `fasta` may be used to store *additional* genome files (not part of in the [RefSeq](https://ftp.ncbi.nlm.nih.gov/refseq/release/)) that should be considered for the database generation process (see [...]()).
* `db` is where Genestrip puts the generated database. If your focus is on filtering fastq files, Genestrip can create specialized, much smaller filtering databases that will be put there too. Moreover, an intermediate database `<project_name>_tempdb.ser.gz` will be put there as part of the database generation process. The intermediate database is not required for *k*-mer matching or fastq filtering and so, it may be deleted afterwards.
* `krakenout` is for output files in the style of [Kraken](https://ccb.jhu.edu/software/kraken/MANUAL.html#output-format). They may optionally be generated when analyzing fastq files.

# Analyzing fastq files by matching *k*-mers of reads

Genestrip's main purpose is to analyze reads from fastq files and count the contained *k*-mers per tax id according to a previously generated database. As an example, Genestrip comes with a small fastq-file `sample.fastq.gz` in `./data/projects/human_virus/fastq`. To start the matching process for it, run
```
sh ./bin/genestrip.sh -f ./data/projects/human_virus/fastq/sample.fastq.gz human_virus match
```
The resulting CSV file will be named `human_virus_match_sample.csv` under `./data/projects/human_virus/csv`.

Here is an example line of it contents along with the header line:
```
name;rank;taxid;reads;kmers;unique kmers;contigs;average contig length;max contig length; max contic desc.;normalized kmers; exp. unique kmers; unique kmers / exp.; quality prediction;
Human gammaherpesvirus 4;NO_RANK;10376;5;158;155;8;19.7500;94;@NS500362:54:HT523BGX2:4:13611:9451:14208 2:N:0:2;0.009677746780291664;157.9111;0.9816;0.00949933507052909;
```
TODO: **Explain meaning of entries here...*

The same principles apply to your own projects under `./data/projects`.

# Filtering fastq files

Genestrip can be used to *filter* fastq files via *k*-mers present in a previously generated database. As an example, one may also use the fastq file `sample.fastq.gz` from `./data/projects/human_virus/fastq`. To start the filtering process, run
````
sh ./bin/genestrip.sh -f ./data/projects/human_virus/fastq/sample.fastq.gz human_virus filter
````
First, the command creates a filtering database file `human_virus_index.ser.gz`, if not yet present, under `./data/projects/<project_name>/db`. 

The resulting filtered fastq file named `human_virus_filtered_sample.fastq.gz` will be put under `./data/projects/human_virus/fastq`. It holds just the reads which contain at least one *k*-mer from the `human_virus` database. Note that, in case of the specific sample fastq file, the filtered fastq file has about the same size as the original file because the sample had *originally been created to just contain human virus-related reads*.

A filtering database is typically even smaller than a database required for *k*-matching but the former can only be used for filtering. So, if you are tight on resources and your focus is on filtering, you may prefer using a filtering database. Also, the filtering process is even faster than the *k*-mer matching process.

# Usage and Goals

The usage of Genestrip:
```
usage: genestrip [options] <project> [<goal1> <goal2>...]
 -d <base dir> Base directory for all data files. The default is
               './data'.
 -f <fqfile>   Input fastq file in case of filtering or k-mer matching
               (regarding goals 'filter' and 'match') or csv file with a
               list of fastq files to be processed via 'multimatch'. Each
               line should have the format '<prefix> <path to fastq
               file>'.
 -r <path>     Common store folder for filtered fastq files and csv files
               created via the goals 'filter' and 'match'. The default is
               '<base dir>/projects/<project name>/fastq' and '<base
               directory>/projects/<project name>/csv' respectively.
 -t <target>   Generation target ('make', 'clean', 'cleanall'). The
               default is 'make'.
```

Genestrip follows a goal-oriented approach in order to create any result file (in similarity to  [make](https://pubs.opengroup.org/onlinepubs/9699919799/utilities/make.html)). Goals are executed in a lazy manner, i.e. a file is only (re-)generated, if it is missing at its designated place in the `<base dir>` folder or any of its subfolders.

Most user-related goals are project-oriented. This means that a database project must exist for the goal to be executable. Some goals refer to fastq files. This means that they can only be executed if a fastq file or a multi-match CSV file is given via the `-f` option. 

Genestrip supports matching of *k*-mers from multiple fastq files in batches: For this purpose you may create a mutli-match CSV file with one line per fastq file to be processed. Each line should have the following form:
```
<name> <path_to_fastq_file>
```
If several fastq files in the multi-match CSV file are associated with the same `<name>`, then the matching results of these files will be aggregated. A resulting CSV file named `<project_name>_multimatch_<name>.csv` will be put under `<base dir>/projects/<project_name>/csv` unless specified otherwise via the `-r` option. If `<path_to_fastq_file>` is a file name without a path prefix, the file is assumed to be located in `<base dir>/projects/<project_name>/fastq`.

Fastq files and fasta file may be g-zipped or not. Genestrip will automatically recognize g-zipped files via the suffixes `.gz` and `.gzip`.

Some named goals are for internal purposes only. In principle they could be run directly by users but rather serve the generation process of databases or they exist for experimental reasons.

Here is the list of user-related goals:
- `show`: Show all goals. Note that some goals will only be shown when using the `-f` option.
- `db`: Generate the database for *k*-mer matching with respect to the given project.
- `dbinfo`: Write information about a project's database content to a CSV file. 
- `index`: Generate a filtering database with respect to a given project.
- `genall`: Generate the *k*-mer matching database and also the filtering database with respect to the given project.
- `clear`: Clear the folders `csv`, `db` and `krakenout`  of a project. This will delete all files the respective folders!
- `match`: Analyze a fastq file as given by the `-f` option. The resulting CSV file will be stored in `<base dir>/projects/<project_name>/csv` unless specified otherwise via the `-r` option.
- `multimatch`: Analyze several fastq files as specified via multi-match CSV file given by the `-f` option.
- `filter`:  Filter a fastq file as given by the `-f` option. The resulting filtered fastq file `filtered_<fqfile>` will be stored under `<base dir>/projects/<project_name>/fastq/` unless specified otherwise via the `-r` option.

Many goals depend on other goals. E.g., the `dbinfo` goal requires the corresponding database to exist and so, it will trigger the execution of the ``db`` goal in case the corresponding database is missing and so on.

# Targets

Genestrip supports three targets for each goal, namely `make`, `clean` and `cleanall`.

- `make` is the default target. It executes a given goal as explained before and creates the files associated with the goal in a lazy manner. If the corresponding files already exist, then `make` does nothing.
- `clean` *deletes* all files associated with the given goal but it does not delete any files of goals that the given goal depends on.
- `cleanall` does same as `clean`, but it *also recursively deletes* any files of goals that the given goal depends on.

# Manually adding fasta-files

In some cases you may want to add *k*-mers of genomes to your database, where the genomes are not part of the [RefSeq](https://ftp.ncbi.nlm.nih.gov/refseq/release/). Genestrip supports this via an optional text file `additional.txt` under `<base dir>/projects/<project_name>`.
The file should contain one line for each additional genome file. (The file must be in fasta format and may be g-zipped or not.)
The line format is
```
<tax id> <path_to_fasta_file>
```
where `<tax id>` is the (unique) tax id associated with the file's genomic data. If `<path_to_fastq_file>` is a file name without a path prefix, then the file is assumed to be located in `<base dir>/projects/<project_name>/fasta`.

# Configuration properties

An optional configuration properties file `Config.properties` may be put under `<base dir>`.
Entries per line should have the form
```
<key>=<value>
```
The following entries are possible:

| Key         | Default Value     | Description |
| ----------- | ----------- | ----------- |
| logLevel      | info       | The levels `error` and `warn` are also possible and result in (much) less verbose logging.       |
| threads      | -1       | The number of consumer threads *n* when processing data with respect to goals ``match``, ``filter`` and ``multimatch`` and also so during the update phase of the ``db`` goal. There is always one additional thread that reads and uncompresses a corresponding fastq or fasta file (so it is *n + 1* threads in total). When negative, the number of available processors *- 1* is used as *n*. When 0, then the corresponding goals run in single-threaded mode. |
| countUniqueKMers      | true       | If `true`, unique kmers will be counted. This requires less than 5% of additional main memory.        |
| writeDumpedFastq   | false        | If `true`, then ``filter`` will also generate a fastq file `dumped_<fqfile>` with all reads not written to the corresponding filtered fastq file. |
| writeFilteredFastq   | false        | If `true`, then the goal `match` writes a filtered fastq file in the same way that the goal `filter` does. Moreover, Genestrip will write an output file `<fqfile>.out` in the [Kraken output format](https://ccb.jhu.edu/software/kraken/MANUAL.html#output-format) under `<base dir>/projects/<project_name>/krakenout` covering all filtered reads. |
| ignoreMissingFastas   | true        | If `true`, then a download of files from NCBI will not stop in case a file is missing on the server.        |
| completeGenomesOnly   | false        | If `true`, then only genomic accessions with the prefixes `AC`, `NC_`, `NZ_` will be considered when generating a database. Otherwise, all genomic accessions will be considered. See [RefSeq accession numbers and molecule types](https://www.ncbi.nlm.nih.gov/books/NBK21091/table/ch18.T.refseq_accession_numbers_and_mole/) for details.       |
| matchWithKMerCounts   | false        | Experimental: Counts how many times each unique *k*-mer has been detected.        |
| maxKMerResCounts   | 200        | The number of the most frequent *k*-mers that will be reported, if `matchWithKMerCounts=true`.       |
| kMerSize   | 31        | The number of base pairs *k* for a *k*-mers. Changes to this values do *not* affect the memory usage of database. A value > 32 will cause collisions, i.e. leads to false positives for the `match` goal. |
| useHttp | true | Use http(s) to download data from NCBI. If ``false``, then Genestrip will try anonymous FTP instead (with login and password set to `anonymous`). 
| refseqHttpBaseURL   | https://ftp.ncbi.nlm.nih.gov/refseq | This [mirror](https://www.funet.fi/pub/mirrors/ftp.ncbi.nlm.nih.gov/refseq/) might be considered as an alternative. (No other mirror sites are known.) |
| refseqFTPBaseURL   | ftp.ncbi.nih.gov       |         |
| taxHttpBaseURL   | https://ftp.ncbi.nlm.nih.gov        | This base URL will be extended by the path `/pub/taxonomy/` in order to download the taxonomy file `taxdmp.zip`.        |
| taxFTPBaseURL   | ftp.ncbi.nih.gov        |         |


# API-based usage

TODO