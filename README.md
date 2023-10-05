Genestrip - Hyper-efficient unique *k*-mer counting and sequence filtering for selected groups of species
===============================================

# Introduction

Metagenomic analysis has become an extremely import field in bio informatics. To analyse large sets of reads, researchers use highly efficient software tools based on *k*-mer matching and counting such as [Kraken](https://github.com/DerrickWood/kraken), [Kraken 2](https://github.com/DerrickWood/kraken2), [KrakenUniq](https://github.com/fbreitwieser/krakenuniq) or [Kaiju](https://github.com/bioinformatics-centre/kaiju). With regard to pathogen detection [KrakenUniq](https://github.com/fbreitwieser/krakenuniq) is particularly suited because of its
 [very low false positive rate](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-018-1568-0) but also because its efficiency, precision and sensitivity.

To avoid false positive classifications of reads or mismatches of *k*-mers, all of these tools resort to very large databases, containing millions of encoded *k*-mers along with their respective tax ids. A database is usually loaded entirely into main memory and consumes tens of gigabytes of space. This requirement mandates specialized and expensive compute servers with very large main memory.

If only a smaller group of a few 100 species or strains is to be considered, one could simply *generate a smaller database that contains just the k-mers of the respective genomes*. **However, when done improperly, this may lead to many false positives at matching time.** The reason for this is that due to the evolutionary relationship of organisms, *a large fraction of k-mers that belong to the genome of a species are unspecific* in a sense that they can as well be found in the genome other species not considered for analysis.

Genestrip offers an efficient and sophisticated database generation process that accounts for this problem: 
* Genstrips's databases can be generated in a fully automated way on the basis of the [RefSeq](https://ftp.ncbi.nlm.nih.gov/refseq/release/). 
* Species or strains, whose *k*-mers should be contained, are specified via a text file containing the corresponding [tax ids](https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/). 
* After the database generation has completed, it can be used for highly efficient *k*-mer matching.
* As a Genestrip database comprises just the genomes of the specified tax ids, it tends to remain small. Therefore analysis can be performed on edge machines as well as on small sized PCs.
* Much like KrakenUniq, Genestrip supports (exact) unique *k*-mer counting and total *k*-mer counting. All match results are stored in a CSV file. 
* Moreover, fastq files can be *filtered* via *k*-mers from a database to greatly reduce their size. Still, all reads containing *k*-mers of respective tax ids will be kept in the filtered fastq file. A similar functionality can be achieved via [biobloom](https://github.com/bcgsc/biobloom), but preparing and applying related filters is very convenient with Genestrip.

# License

[Genestrip is free for non-commercial use.](./LICENSE.txt) Please contact [daniel.pfeifer@progotec.de](mailto:daniel.pfeifer@progotec.de) if you are interested in a commercial license.

# Building and installing

Genestrip is structured as a standard [Maven 2 or 3](https://maven.apache.org/) project and is compatible with the [JDK](https://jdk.java.net/) 1.8 or higher.[^1]

`cd` to the installation directory `genestrip`. Given a matching Maven and JDK installation,
`mvn install` will compile, test and install the Genestrip program library.
Afterwards a self-contained and executable `genestrip.jar` file will be stored under `./lib`. 

[^1]: Counter to common belief Java, can well be used for such high performance applications when using its programming facilities the right way.

# Building the sample database

The Genestrip installation holds additional folders that include the sample project `human_virus`.
After building Genestrip, you may call
`./bin/genestrip.sh human_virus storeinfo`
in order to make the `human_virus` database and create CSV file with basic information about the database content.

Genestrip follows a goal oriented approach in order to create any result files (in the spirit of  [make](https://pubs.opengroup.org/onlinepubs/9699919799/utilities/make.html)). So, when building the `human_virus` database, Genestrip will:
1. download the [taxonomy](https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip) and unzip it to `./data/common`,
1. download the [refseq release catalog](https://ftp.ncbi.nlm.nih.gov/refseq/release/release-catalog/) to `./data/common/refseq`,
1. download [all virus related RefSeq fna files](https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/) to `./data/commmon/refseq` (which is currently just a single file),
1. perform several follow-up goals until the database file `human_virus_store.ser.gz` is finally made under `./data/projects/human_virus/filters`, 
1. create a CSV file `human_virus_storeinfo.csv` under `./data/projects/human_virus/csv` which contains basic information about the database, i.e. the number of *k*-mers stored per tax id.

The generated database comprises *k*-mers for all viruses according to the tax id file `./data/project/human_virus/taxids.txt`.


# Building your own database

Building your own database is straight-forward:
1. Create a folder `<db_name>` under `./data/projects/`. The folder is the places for all subsequent files to specify the content of the database. It is also the core name of the database.
1. Create a text file `./data/projects/<db_name>/taxids.txt` with one tax id per line. The database will *only* contain *k*-mers from genomes that belong to tax ids from the file `taxids.txt`. If `taxids.txt` contains a tax id from a non-leaf node of the [taxonomy](https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip), then all subordinate tax ids and *k*-mers from respective [RefSeq](https://ftp.ncbi.nlm.nih.gov/refseq/release/) genomes will be included in the database.
1. Create a text file `./data/projects/<db_name>/categories.txt`. The file tells Genestrip, which categories of organisms should be *considered* when building the database and when determining a least common ancestor tax id for *k*-mer in the finalized database. You also have to ensure that the categories from `categories.txt` comprise all of the tax ids from `taxids.txt`: Only genomes of comprised tax ids will be used to generate the database.
1. The following categories are allowed, and originate from the [RefSeq](https://ftp.ncbi.nlm.nih.gov/refseq/release/):
`archaea`, `bacteria`, `complete`, `fungi`, `invertebrate`, `mitochondrion`, `other`, `plant`, `plasmid`, `plastid`, `protozoa`, `vertebrate_mammalian`, `vertebrate_other` and `viral`. You may enter one category per line.
The *entire* set of genomes that belong a category referenced in `categories.txt` will be used to determine the least common ancestors tax ids of the *k*-mers stored in the finalized database. The more categories you reference in `categories.txt`, the more genomes files will have to be downloaded from the [RefSeq](https://ftp.ncbi.nlm.nih.gov/refseq/release/) and the longer the database generation process will take. E.g., if only `viral` is included, it will take only a few minutes, but with `bacteria` included, it may typically take about a day or more.
So you should chose a suitable set of categories to balance the completeness of considered genomes and efficiency of the database generation process.

1. Start the database generation process via the command `./bin/genestrip.sh <db_name> storeinfo`. This will also create a CSV file `./data/projects/<db_name>/csv/<db_name>_storeinfo.csv` with basic information about the database, i.e. the number of *k*-mers stored per tax id. The path of the database will be `./data/projects/<db_name>/filters/<db_name>_store.ser.gz`

In general, Genestrip organizes a project folder `./data/projects/<db_name>` by means of the following sub-folders:
* `csv` is where analysis results of fastq files will be stored (by default).
* `fastqs` is where Genestrip will store filtered fastq files. You may also place the fastq files to be analyzed there (but they can be read from any other path too).
* `fastas` can be used to store *additional* genome files (not part of in the [RefSeq](https://ftp.ncbi.nlm.nih.gov/refseq/release/)) that should be considered for database generation process (see [...]()).
* `filters` is where Genestrip places the generated database. If your focus is on filtering fastq files, Genestrip can create specialized, much smaller filtering databases that will be placed there too. Moreover an interim database `<db_name>_tempstore.ser.gz` will be put there as part of the database generation process. It is not required for *k*-mer matching or fastq filtering and so, it may be deleted afterwards.
* `krakenout` is for output files in the style of [Kraken](https://ccb.jhu.edu/software/kraken/MANUAL.html#output-format). They may optionally be generated when analyzing fastq files.

# Matching *k*-mers of reads from fastq-files

Based on a database Genestrip can analyze reads from fastq-file and count the contained *k*-mers according their tax ids.
Genestrip comes with a smalle sample fastq-file `demo.fastq.gz` in `./data/projects/human_virus/fastq`. To start the matching process for it, run
```
./bin/genestrip.sh -f ./data/projects/human_virus/fastq/demo.fastq.gz human_virus match
```
The resulting CSV file will be named `human_virus.demo.fastq.gz.csv`, is stored in `./data/projects/human_virus/csv` and should have the following content:
```
```

# Goals

# Configuration details

# Manually adding fasta-files

# API-based use

