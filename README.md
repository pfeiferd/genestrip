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

# Building a sample database

# Matching *k*-mers from fastq-files

# Building you own database

# Goals

# Configuration details

# Manually adding fasta-files

# API-base use 