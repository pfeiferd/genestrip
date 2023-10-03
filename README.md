Genestrip - Hyper-efficient unique *k*-mer counting and sequence filtering for selected groups of species
===============================================

# Introduction

Metagenomic analysis has become an extremely import field in bio informatics. To analyse large sets of reads researchers use highly efficient software tools based on *k*-mer matching and counting such as [Kraken](https://github.com/DerrickWood/kraken), [Kraken 2](https://github.com/DerrickWood/kraken2), [KrakenUniq](https://github.com/fbreitwieser/krakenuniq) or [Kaiju](https://github.com/bioinformatics-centre/kaiju). With regard to pathogen detection [KrakenUniq](https://github.com/fbreitwieser/krakenuniq) is particularly suited because of its
 [very low false positive rate]() but also because its efficiency, precision and sensitivity.

To avoid false positive classifications of reads or mismatches of *k*-mers, all of these tools resort to very large databases, containing millions of encoded *k*-mers along with their respective tax ids. A database is usually loaded entirely into main memory and consumes tens of gigabytes of space. This requirement mandates specialized and expensive compute servers with very large main memory.

If only a smaller group of a few 100 species or strains is to be considered, one may simply *generate a smaller database that contains just the k-mers of the respective genomes*. **However, when done improperly, this may lead to many false positives at matching time.** The reason for this is that due to the evolutionary relationship of organisms, *a large fraction of k-mers that belong to the genome of a species are unspecific* in a sense that they can as well be found in the genome other species not considered for analysis.

Genestrip offers an efficient and sophisticated database generation process that accounts for this problem: 
* Genstrips's databases can be generated in a fully automated way on the basis of the [RefSeq](https://ftp.ncbi.nlm.nih.gov/refseq/release/). 
* Species or strains, whose *k*-mers should be contained, are specified via a text file containing the corresponding [tax ids](https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/). 
* After the database generation has completed, it can be used for highly efficient *k*-mer matching.
* Much like KrakenUniq, Genestrip supports (exact) unique *k*-mer counting and total *k*-mer counting. All match results are stored in a CSV file. 
* Moreover fastq files can be filtered via *k*-mers from a database to greatly reduce their size but by keeping all reads containing respective *k*-mers. This functionality is similar to [biobloom], but preparing a related filter much more sophisticated in Genestrip.

# License

[Genestrip is free for non-commercial use.](./LICENSE.txt) Please contact [daniel.pfeifer@progotec.de](mailto:daniel.pfeifer@progotec.de) if you are interested in a commercial license. Note that illegal use of this software will be prosecuted.

# Building and installing

Genestrip is structured as a standard [Maven 2 or 3](https://maven.apache.org/) project and is compatible with [JDK](https://jdk.java.net/) 1.8 or higher.

Given a Maven and JDK installation,
`mvn install` will compile, test and install the genestrip program library.
Afterwards a self-contained and executable `genestrip.jar` file will stored under `./lib`. 

# Building a sample database

Genestrip comes with a folder structure which includes the project `human_virus`.
After building genestrip, you may call
`./bin/genestip.sh -d ./data human_virus storeinfo`
in order to build the `human_virus` database `human_virus.ser.gz` under `./data/project/human_virus/filters`

Genestrip follows a goal oriented approach in order to create any result files much like [make](https://pubs.opengroup.org/onlinepubs/9699919799/utilities/make.html). So to build the human virus database, genestrip will
1. Download the [taxonomy] and unzip it to `./data/common.`
1. Download the [refseq release catalog[ to `./data/refseq`
1. Download [all virus related RefSeq fna files] to `./data/commmon/refseq`
1. Perform several follow up goals until the `human_virus` database is finally made. 
1. Create a CSV file `storeinfo.csv` under `./data/human_virus/csv` which contains basic information about the database, i.e. the number of $k$-mers stored per taxid.

The database comprises $k$-mers for all viruses according to the tax id file `./data/project/human_virus/taxids.txt`.


# Matching *k*-mers from fastq-files

# Building you own database

# Goals

# Configuration details

# Manually adding fasta-files

# API-base use 