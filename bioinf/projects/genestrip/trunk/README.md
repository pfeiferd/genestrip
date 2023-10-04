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
* As a genestrip database comprises just the genomes of the specified tax ids it tends to remain small. Therefore analysis can be performed in edge machines or on small sized PCs.
* Much like KrakenUniq, Genestrip supports (exact) unique *k*-mer counting and total *k*-mer counting. All match results are stored in a CSV file. 
* Moreover, fastq files can be filtered via *k*-mers from a database to greatly reduce their size but by keeping all reads containing *k*-mers of respective tax ids. A similar functionality can be achieved via [biobloom](https://github.com/bcgsc/biobloom), but preparing related filters is more convenient with Genestrip.

# License

[Genestrip is free for non-commercial use.](./LICENSE.txt) Please contact [daniel.pfeifer@progotec.de](mailto:daniel.pfeifer@progotec.de) if you are interested in a commercial license.

# Building and installing

Genestrip is structured as a standard [Maven 2 or 3](https://maven.apache.org/) project and is compatible with the [JDK](https://jdk.java.net/) 1.8 or higher.[^1]

Given a Maven and JDK installation,
`mvn install` will compile, test and install the genestrip program library.
Afterwards a self-contained and executable `genestrip.jar` file will stored under `./lib`. 

[^1]: Counter to common belief Java can be used for such high performance applications when using it's programming facilities the right way.

# Building a sample database

Genestrip comes with a folder structure which includes the project `human_virus`.
After building genestrip, you may call
`./bin/genestip.sh -d ./data human_virus storeinfo`
in order to build the `human_virus` database `human_virus.ser.gz` under `./data/project/human_virus/filters`

Genestrip follows a goal oriented approach in order to create any result files much like [make](https://pubs.opengroup.org/onlinepubs/9699919799/utilities/make.html). So to build the human virus database, genestrip will:
1. download the [taxonomy](https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip) and unzip it to `./data/common.`,
1. download the [refseq release catalog](https://ftp.ncbi.nlm.nih.gov/refseq/release/release-catalog/) to `./data/refseq`,
1. download [all virus related RefSeq fna files](https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/) to `./data/commmon/refseq` (which is currently just one),
1. perform several follow up goals until the database `human_virus.ser.gz` database is finally made under `./data/projects/human_virus/filters`, 
1. create a CSV file `storeinfo.csv` under `./data/human_virus/csv` which contains basic information about the database, i.e. the number of *k*-mers stored per tax id.

The generated database comprises *k*-mers for all viruses according to the tax id file `./data/project/human_virus/taxids.txt`.


# Matching *k*-mers from fastq-files

# Building you own database

# Goals

# Configuration details

# Manually adding fasta-files

# API-base use

