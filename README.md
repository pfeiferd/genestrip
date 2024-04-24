
[comment]: # (“Commons Clause” License Condition v1.0)
[comment]: # ()  
[comment]: # (The Software is provided to you by the Licensor under the License,)
[comment]: # (as defined below, subject to the following condition.)
[comment]: # ()  
[comment]: # (Without limiting other conditions in the License, the grant of rights under the License)
[comment]: # (will not include, and the License does not grant to you, the right to Sell the Software.)
[comment]: # ()
[comment]: # (For purposes of the foregoing, “Sell” means practicing any or all of the rights granted)
[comment]: # (to you under the License to provide to third parties, for a fee or other consideration)
[comment]: # (including without limitation fees for hosting or consulting/ support services related to)
[comment]: # (the Software, a product or service whose value derives, entirely or substantially, from the)
[comment]: # (functionality of the Software. Any license notice or attribution required by the License)
[comment]: # (must also include this Commons Clause License Condition notice.)
[comment]: # ()
[comment]: # (Software: genestrip)
[comment]: # ()
[comment]: # (License: Apache 2.0)
[comment]: # ()
[comment]: # (Licensor: Daniel Pfeifer, daniel.pfeifer@progotec.de)

Genestrip - Efficient read classification, filtering and *k*-mer counting for selected groups of species 
===============================================

# Introduction

Metagenomic analysis has become an extremely import field in bio informatics. To analyse large sets of reads, researchers use highly efficient software tools based on *k*-mer matching and counting such as [Kraken](https://github.com/DerrickWood/kraken), [Kraken 2](https://github.com/DerrickWood/kraken2) or [KrakenUniq](https://github.com/fbreitwieser/krakenuniq). With regard to pathogen detection [KrakenUniq](https://github.com/fbreitwieser/krakenuniq) is particularly suited because of its
 [very low false positive rate](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-018-1568-0) but also because its efficiency, precision and sensitivity.

To avoid false positive classifications of reads or mismatches of *k*-mers, all of these tools resort to very large databases, containing millions of encoded *k*-mers along with their respective tax ids. A database is usually loaded entirely into main memory and consumes tens of gigabytes of space. This requirement mandates specialized and expensive compute servers with very large main memory.

If only a small group of a few dozen species or strains is to be considered, one could simply *generate a smaller database that contains just the *k*-mers of the respective genomes*. **However, when done improperly, this may lead to many false positives at read analysis time.** The reason for this is that due to the evolutionary relationship of organisms, *a large fraction of k-mers that belong to the genome of a species are unspecific* in a sense that they can as well be found in the genome other species not considered for analysis.

Genestrip offers an efficient and sophisticated database generation process that accounts for this problem: 
* Genstrips's databases can be generated in a fully automated way on the basis of the [RefSeq](https://ftp.ncbi.nlm.nih.gov/refseq/release/). On top, genomic files from [Genbank](https://ftp.ncbi.nlm.nih.gov/genomes/genbank/) can be automatically included for the selected group of species or strains. Including such files may refine the database by adding species-related *k*-mers not found in the RefSeq alone.
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
1. download the [RefSeq release catalog](https://ftp.ncbi.nlm.nih.gov/refseq/release/release-catalog/) to `./data/common/refseq`,
1. download [all virus related RefSeq fna files](https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/) to `./data/commmon/refseq` (which is currently just a single file),
1. perform several follow-up goals, until the database file `human_virus_db.kmers.ser.gz` is finally made under `./data/projects/human_virus/db`, 
1. create a CSV file `human_virus_dbinfo.csv` under `./data/projects/human_virus/csv`, which contains basic information about the database, i.e. the number of *k*-mers stored per tax id.

The generated database comprises *k*-mers for all viruses according to the tax id file `./data/project/human_virus/taxids.txt`.


# Generating your own database

Generating your own database is straight-forward:
1. Create a project folder `<project_name>` under `./data/projects/`. This is the place for all subsequent files that specify the content of a database to be generated. It is also the core name of a related database.
1. Create a text file `./data/projects/<project_name>/taxids.txt` with one tax id per line. The database will *only* contain *k*-mers from genomes that belong to tax ids referenced in the file `taxids.txt`. If `taxids.txt` contains a tax id from a non-leaf node of the [taxonomy](https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip), then all subordinate tax ids and *k*-mers from respective [RefSeq](https://ftp.ncbi.nlm.nih.gov/refseq/release/) genomes will be included in the database as well.
1. Create a text file `./data/projects/<project_name>/categories.txt`. This file tells Genestrip, which categories of organisms should be *considered* when generating the database and when determining a least common ancestor tax id for *k*-mers in the finalized database. You also have to ensure that the categories from `categories.txt` comprise all of the tax ids from `taxids.txt`: Only genomes of comprised tax ids will be used to generate the database.
1. The following categories are allowed, and originate from the [RefSeq](https://ftp.ncbi.nlm.nih.gov/refseq/release/):
`archaea`, `bacteria`, `complete`, `fungi`, `invertebrate`, `mitochondrion`, `other`, `plant`, `plasmid`, `plastid`, `protozoa`, `vertebrate_mammalian`, `vertebrate_other` and `viral`. You may enter one category per line.
The *entire* set of genomes that belong to a category referenced in `categories.txt` will be used to determine the least common ancestors tax ids of the *k*-mers stored in the finalized database. The more categories you reference in `categories.txt`, the more genomes files will have to be downloaded from the [RefSeq](https://ftp.ncbi.nlm.nih.gov/refseq/release/) and the longer the database generation process will take. E.g., if only the category `viral` is included, it should just take a few minutes, but with `bacteria` included, it may typically take about a day or more.
So you should chose a suitable set of categories to balance completeness of considered genomes with efficiency of the database generation process.

1. Start the database generation process via the command `sh ./bin/genestrip.sh <project_name> dbinfo`. This will also create a CSV file `./data/projects/<project_name>/csv/<project_name>_dbinfo.csv` with basic information about the database, i.e. the number of *k*-mers stored per tax id. The path of the database will be `./data/projects/<project_name>/db/<project_name>_db.kmers.ser.gz`. The file `./data/projects/<project_name>/db/<project_name>_db.bloom.ser.gz` is an optional part of the database and created as well. It is used by default during fastq file analysis to improve performance (in most cases) but increases main memory usage. (Please check the configuration property `useBloomFilterForMatch` for more information on this.)

In general, Genestrip organizes a project folder `./data/projects/<project_name>` by means of the following sub-folders:
* `csv` is where analysis results of fastq files will be stored (by default).
* `fastq` is where Genestrip will store filtered fastq files. You may also put the fastq files to be analyzed there (but they can be read from any other path too).
* `fasta` may be used to store *additional* genome files (not part of in the [RefSeq](https://ftp.ncbi.nlm.nih.gov/refseq/release/)) that should be considered for the database generation process (see Section [Manually adding fasta-files](#manually-adding-fasta-files)).
* `db` is where Genestrip puts the generated database. If your focus is on filtering fastq files, Genestrip can create a specialized, smaller filtering database named `<project_name>_db.bloom.ser.gz` that will be put there too. Moreover, the intermediate databases `<project_name>_tempdb.kmers.ser.gz` and `<project_name>_tempdb.bloom.ser.gz` will be put there as part of the database generation process. The intermediate databases are not required for *k*-mer matching or fastq filtering and so, they may be deleted afterwards.
* `krakenout` is for output files in the style of [Kraken](https://ccb.jhu.edu/software/kraken/MANUAL.html#output-format). They may optionally be generated when analyzing fastq files.

A database covering more species may require more memory - especially when generating the database. This can be addressed by adjusting the script `genestrip.sh` where the argument `-Xmx32g` sets the maximum heap space of the Java Virtual Machine to 32GB. E.g. to double it, simple replace `32g` by `64g`.


# Analyzing fastq files by matching *k*-mers of reads

Genestrip's main purpose is to analyze reads from fastq files and count the contained *k*-mers per tax id according to a previously generated database. As an example, Genestrip comes with a small fastq-file `sample.fastq.gz` in `./data/projects/human_virus/fastq`. To start the matching process for it, run
```
sh ./bin/genestrip.sh -f ./data/projects/human_virus/fastq/sample.fastq.gz human_virus match
```
The resulting CSV file will be named `human_virus_match_sample.csv` under `./data/projects/human_virus/csv`. The same principles apply to your own projects under `./data/projects`.

These are a few lines of its contents along with the header line:
```
name;rank;taxid;reads;kmers from reads;kmers;unique kmers;contigs;average contig length;max contig length;max contig desc.;normalized kmers;exp. unique kmers;unique kmers / exp.;quality prediction;
TOTAL;SUPERKINGDOM;1;6565;0;461305;0;0;0;0;0;0;0;0;0;
Viruses;SUPERKINGDOM;10239;1016;15939;27712;72;2653;10.4568;101;@NS500362:54:HT523BGX2:4:12412:20538:2021 2:N:0:2;2960.144020983948;80.0000;0.9000;2664.129618885553;
Orthornavirae;KINGDOM;2732396;19;159;122;4;119;1.2773;32;@NS500362:54:HT523BGX2:4:11405:19113:11578 2:N:0:2;260.63623741342496;4.0000;1.0000;260.63623741342514;
Punta Toro virus;NO_RANK;11587;27;50;50;12;30;2.6667;36;@NS500362:54:HT523BGX2:4:11502:8048:10170 2:N:0:2;0.03412447251416636;49.9023;0.2405;0.008205909512196132;
Ross River virus;SPECIES;11029;121;889;195;7;124;1.8145;35;@NS500362:54:HT523BGX2:4:23501:24909:11663 2:N:0:2;0.08970031917772085;193.9853;0.0361;0.0032368543472774043;
```
The meaning of the columns is a follows:

| Column      | Description |
| ----------- | ----------- |
| `name`      | The name associated with the tax id. |
| `rank`      | The rank of the tax id. |
| `taxid`      | The tax id. |
| `reads`      | The mumber of reads classified with respect to the tax id. |
| `kmers from reads`      | The number of *k*-mers from classified reads which are consistent with the read's tax id. |
| `kmers`      | *All* matched *k*-mers which are specific to the tax id's genome (according to the database). The *k*-mers do not have to be in an accordingly classified read for this. |
| `unique kmers` | *All* unique *k*-mers, which are specific to the tax id's genome (according to the database). Here, multiple occurrences of the same *k*-mer are only counted once. Genestrip always performs exact counting according to [KrakenUniq's exact counting](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-018-1568-0#Sec8). (Genestrip implements an efficient in-memory storage method for related counts based on a bit vectors.) |
| `contigs`      |  The number of contiguous sequences of *k*-mers that are specific to the tax id's genome.     |
| `average contig length`      | The average  base pair length of contiguous sequences of *k*-mers that are specific to the tax id's genome. |
| `max contig length`      |  The maximum base pair length of all contiguous sequences of *k*-mers that are specific to the tax id's genome. |
| `max contig desc.`      |  The descriptor of a read that holds a contiguous sequence of maximum length (according to the previous column).   |
| `normalized kmers`      |   *k*-mer counts from column `kmers` but normalized with respect to the total number of *k*-mers per fastq file and the number of specific *k*-mers for a tax id in the database. The value allows for a less biased comparison of *k*-mer counts across fastq files and across species. It is computed as `normalizedKMersFactor` * `kmers` */ k<sub>f</sub> * u / u<sub>t</sub>*, where *k<sub>f</sub>* is the total number of *k*-mers in the fastq file, *u* is the total number of *k*-mers in the database and *u<sub>t</sub>* is the number of specific *k*-mers for a tax id in the database. `normalizedKMersFactor` is a configuration property; its default is 1000000000 (see also Section [Configuration properties](#configuration-properties)). |
| `exp. unique kmers`      |  The number of expected unique *k*-mers, which is *u<sub>t</sub> * (1 - (1 - 1/u<sub>t</sub>)*<sup>`kmers`</sup>), where *u<sub>t</sub>* is the number of specific *k*-mers for a tax id in the database. |
| `unique kmers / exp.`      |  The ratio `unique kmers` / `exp. unique kmers` for the tax id. This should be close to 1 for a consistent match of *k*-mers. ([This paper](https://arxiv.org/pdf/1602.05822.pdf) discusses the corresponding background distribution (of `unique kmers`).)  |
| `quality prediction`      | Computed as  `normalized kmers` * `unique kmers / exp.`. It combines the normalized counts of *k*-mers with the valued consistency between *k*-mers and unique *k*-mers. |
| `max kmer counts` | The frequencies of the most frequent unique *k*-mers which are specific to the tax id's genome in descending order separated by `;`. This column is experimental and only added when the configuration property `matchWithKMerCounts` is set to `true`. The number of frequencies is determined via `maxKMerResCounts` (see also Section [Configuration properties](#configuration-properties)). |

The frequencies from `max kmer counts` can be used to build frequency graph's for *k*-mers as shown below. The frequency graphs help to further assess the validity of analysis results.

<p align="center">
  <img src="EBV-Sample-Graph.png" width="800"/>
</p>
The depicted graph is based on the following CSV result file entry:

```
name;rank;taxid;reads;kmers from reads;kmers;unique kmers;contigs;average contig length;max contig length;max contig desc.;db coverage;normalized kmers;exp. unique kmers;unique kmers / exp.;quality prediction;max kmer counts;
...
Human gammaherpesvirus 4;SPECIES;10376;113252;4061786;4151610;120419;239293;47.3495;101;@A01245:81:HTH3LDSX2:4:1248:18231:31093 1:N:0:AATATTGCCA+GGTGTCACCG;0.8631;3769564.5570;139514.0000;0.8631;3253631.8534;629;625;590;453;449;445;442;433;413;412;411;411;409;409;408;407;400;399;399;397;396;396;396;395;394;393;393;390;389;387;387;387;385;384;384;384;384;383;383;383;383;381;381;380;379;379;379;378;378;377;377;377;376;376;376;375;374;374;372;372;371;370;370;369;368;368;368;368;368;367;367;367;367;366;366;365;363;363;362;362;361;360;359;359;358;358;357;356;355;355;354;354;354;354;354;353;353;352;352;352;352;351;351;351;351;350;350;349;349;348;348;348;348;348;347;347;347;347;347;347;345;345;345;345;345;344;344;344;344;344;343;343;343;343;343;343;342;342;342;342;342;341;341;341;341;341;341;340;340;340;340;340;339;339;339;338;338;338;337;337;337;335;335;335;335;335;335;335;334;334;334;334;334;334;333;333;332;332;332;332;331;331;331;331;331;330;330;330;330;330;329;329;329;328;328;328;327;327;327;327;327;327;327;327;327;327;326;326;326;326;326;325;325;325;325;325;325;325;324;324;324;324;324;324;323;323;323;323;323;323;323;323;322;322;322;322;321;321;320;320;320;320;319;319;319;319;319;319;319;319;319;318;318;318;318;317;317;317;317;317;317;317;316;316;316;316;316;316;315;315;315;315;314;314;314;314;314;314;313;313;313;313;313;312;312;312;312;312;312;312;311;311;310;310;310;310;310;310;310;309;309;309;309;309;308;308;308;308;308;308;308;307;307;307;307;307;307;307;307;307;307;306;306;306;306;306;306;306;306;306;306;305;305;305;305;305;305;304;304;304;304;304;303;303;303;303;303;303;303;302;302;302;302;302;302;302;302;302;302;301;301;301;301;301;301;301;301;300;300;300;300;300;300;300;300;300;300;299;299;299;299;299;299;299;299;299;298;298;298;298;298;298;297;297;297;297;297;297;297;297;297;297;297;297;296;296;296;296;296;296;296;296;296;295;295;295;295;295;295;295;295;295;295;294;294;294;294;294;294;294;294;294;294;294;293;293;293;293;293;293;292;292;292;292;292;292;292;291;291;291;290;290;290;290;290;290;289;289;289;289;289;289;289;289;289;289;288;288;288;288;288;288;288;288;288;288;287;287;287;287;287;287;287;287;287;287;287;287;287;286;286;286;286;286;286;286;286;286;286;286;
```
As this result is rather consistent with the statistical expectation (for the unique *k*-mer frequency distribution), the graph is quite flat and `unique kmers / exp. = 0.8631` is close to 1.

# Reliability of results

We cannot guarantee for any results returned by Genestrip. Use this software at you own risk. **Important: It is by no means meant to be used for any medical purposes** and it is purely experimental in nature.

Despite of the these limitations, we tested the Genestrip in the following ways:
* Critical code is covered by functional tests (using [JUnit](https://junit.org)).
* Results have been evaluated based on [Kraken's accuracy datasets](https://ccb.jhu.edu/software/kraken/dl/accuracy.tgz) by including "HiSeq" and "MiSeq". For this purpose we generated Genestrip databases covering the organisms from respective datasets. Genestrip's accuracy results were comparable to those obtained via [Kraken](https://ccb.jhu.edu/software/kraken/).
* We applied the human virus database from above to real-world fastq files, and Genestrip returned similar results and (unique) *k*-mer counts as [KrakenUniq](https://github.com/fbreitwieser/krakenuniq) for this case.
* We applied Genestrip to 20 real-world fastq files based on human saliva samples. The findings matched the "general expectations" with regard to Herpes viruses and mouth bacteria such as Steptococcus mutans, Helicobacter pylori and others.

# Filtering fastq files

Genestrip can be used to *filter* fastq files via *k*-mers present in a previously generated database. As an example, one may also use the fastq file `sample.fastq.gz` from `./data/projects/human_virus/fastq`. To start the filtering process, run
````
sh ./bin/genestrip.sh -f ./data/projects/human_virus/fastq/sample.fastq.gz human_virus filter
````
First, the command creates a filtering database file `human_virus_index.bloom.ser.gz`, if not yet present, under `./data/projects/<project_name>/db`. 

The resulting filtered fastq file named `human_virus_filtered_sample.fastq.gz` will be put under `./data/projects/human_virus/fastq`. It holds just the reads which contain at least one *k*-mer from the `human_virus` database. Note that, in case of the specific sample fastq file, the filtered fastq file has about the same size as the original file because the sample had *originally been created to just contain human virus-related reads*.

A filtering database is typically smaller than a database required for *k*-mer matching but the former can only be used for filtering. So, if you are tight on resources and your focus is on filtering, you may prefer using a filtering database. Also, the filtering process is  faster than the *k*-mer matching process.

# Usage and Goals

The usage of Genestrip:
```
usage: genestrip.sh [options] <project> [<goal1> <goal2>...]
 -d <base dir>   Base directory for all data files. The default is
                 './data'.
 -f <fqfile>     Input fastq file in case of filtering or k-mer matching
                 (regarding goals 'filter' and 'match') or CSV file with a
                 list of fastq files to be processed via 'multimatch'.
                 Each line of a CSV file should have the format '<prefix>
                 <path to fastq file>'.
 -r <path>       Common store folder for filtered fastq files and
                 resulting CSV files created via the goals 'filter' and
                 'match'. The defaults are '<base dir>/projects/<project
                 name>/fastq' and '<base dir>/projects/<project
                 name>/csv', respectively.
 -t <target>     Generation target ('make', 'clean' or 'cleanall'). The
                 default is 'make'.
 -tx <taxids>    List of tax ids separated by ',' (but no blanks) for the
                 goal 'db2fastq'. A tax id may have the suffix '+', which
                 means that taxonomic descendants from a database will be included.
```

Genestrip follows a goal-oriented approach in order to create any result file (in similarity to  [make](https://pubs.opengroup.org/onlinepubs/9699919799/utilities/make.html)). Goals are executed in a lazy manner, i.e. a file is only (re-)generated, if it is missing at its designated place in the `<base dir>` folder or any of its subfolders.

Most user-related goals are project-oriented. This means that a database project must exist for the goal to be executable. Some goals refer to fastq files. This means that they can only be executed if a fastq file or a multi-match CSV file is given via the `-f` option. 

Genestrip supports matching of *k*-mers from multiple fastq files in batches: For this purpose you may create a mutli-match CSV file with one line per fastq file to be processed. Each line should have the following form:
```
<name> <path_to_fastq_file>
```
If several fastq files in the multi-match CSV file are associated with the same `<name>`, then the matching results of these files will be aggregated. A resulting CSV file named `<project_name>_multimatch_<name>.csv` will be put under `<base dir>/projects/<project_name>/csv` unless specified otherwise via the `-r` option. If `<path_to_fastq_file>` is a file name without a path prefix, the file is assumed to be located in `<base dir>/projects/<project_name>/fastq`.

Fastq files and fasta file may be g-zipped or not. Genestrip will automatically recognize g-zipped files via the suffixes `.gz` and `.gzip`.

Some named goals are for internal purposes only. In principle, they could be run directly by users but rather serve the generation process of databases or they exist for experimental reasons.

Here is the list of user-related goals:
- `show`: Show user-related goals. Note that some goals will only be shown when using the `-f` option.
- `showall`: Show user-related and most internal goals. Note that some goals will only be shown when using the `-f` option.
- `db`: Generate the database for *k*-mer matching with respect to the given project.
- `dbinfo`: Write information about a project's database content to a CSV file.
- `db2fastq`: Generate fastq files from the database. A respective fastq file will contain all *k*-mers specifically associated with a 
single tax id from the database where each *k*-mer is represented by a read consisting of *k* bases. Respective fastq files will be stored 
in `<base dir>/projects/<project_name>/fastq` with the file name format `<project_name>_db2fastq_<taxid>.fastq.gz`. 
The command line option `tx` serves at selecting the corresponding tax ids for the fastq files to be generated (see Section [Usage and Goals](#usage-and-goals)). 
If the option is omitted, then fastq files for *all* tax ids from the database will be generated.
- `index`: Generate a filtering database with respect to a given project.
- `genall`: Generate the *k*-mer matching database and also the filtering database with respect to the given project.
- `clear`: Clear the folders `csv`, `db` and `krakenout`  of a project. This will delete all files the respective folders!
- `match`: Analyze a fastq file as given by the `-f` option. The resulting CSV file will be stored in `<base dir>/projects/<project_name>/csv` unless specified otherwise via the `-r` option.
- `multimatch`: Analyze several fastq files as specified via multi-match CSV file given by the `-f` option.
- `matchlr`: Same as `match` but without doing read classification. This corresponds to the configuration setting `classifyReads=false` from below.
- `multimatchlr`:  Same as `multimatch` but without doing read classification. This corresponds to the configuration setting `classifyReads=false` from below.
- `filter`:  Filter a fastq file as given by the `-f` option. The resulting filtered fastq file `filtered_<fqfile>` will be stored under `<base dir>/projects/<project_name>/fastq/` unless specified otherwise via the `-r` option.

Many goals depend on other goals. E.g., the `dbinfo` goal requires the corresponding database to exist and so, it will trigger the execution of the ``db`` goal in case the corresponding database is missing and so on.

# Targets

Genestrip supports three targets for each goal, namely `make`, `clean` and `cleanall`.

- `make` is the default target. It executes a given goal as explained before and creates the files associated with the goal in a lazy manner. If the corresponding files already exist, then `make` does nothing.
- `clean` *deletes* all files associated with the given goal but it does not delete any files of goals that the given goal depends on.
- `cleanall` does same as `clean`, but it *also recursively deletes* any files of goals that the given goal depends on.

# Manually adding fasta files

In some cases you may want to add *k*-mers of genomes to your database, where the genomes are not part of the [RefSeq](https://ftp.ncbi.nlm.nih.gov/refseq/release/). Genestrip supports this via an optional text file `additional.txt` under `<base dir>/projects/<project_name>`.
The file should contain one line for each additional genome file. (The genome file must be in fasta format and may be g-zipped or not.)
The line format is
```
<taxid> <path_to_fasta_file>
```
where `<taxid>` is the (unique) tax id associated with the file's genomic data. (Multiple tax ids per fasta file are not supported in this context.) 
If `<path_to_fastq_file>` is a file name without a path prefix, then the file is assumed to be located in `<base dir>/projects/<project_name>/fasta`. If not found there,
the directory `<base dir>/fasta` will be checked as a secondary location.

This adding of fasta files can also be used to *just* correct the least common ancestor of *k*-mers in the resulting database since the added fasta files will be automatically used during the update phase of the ``db`` goal. E.g., to correct the least common ancestor of *k*-mers occurring in a purely `protozoa`n database *but also* in the human genome, one may simple add
```
9606 <path_to_human_genome_fasta_file>
```
to `additional.txt` where `<path_to_human_genome_fasta_file>` is the path to the [human genome fasta file](https://www.ncbi.nlm.nih.gov/datasets/taxonomy/9606/). (Note that for this purpose, ``9606`` *must not* occur in `taxids.txt`, since otherwise all *k*-mers from the human genome would be included in a `protozoa`n database.)

# Configuration properties

An optional configuration properties file `config.properties` or `Config.properties` may be put under `<base dir>`.
Entries per line should have the form
```
<key>=<value>
```
The following entries are possible:

| Key         | Default Value     | Description | Affected Goals |
| ----------- | ----------- | ----------- | ----------- |
| `logLevel`      | `trace`       | The log levels `error`, `warn` and `info` are also possible and result in (much) less verbose logging. | all |
| `matchLogUpdateCycle`      | `1000000`       | Affects the log level `trace`: Defines after how many reads per fastq file, information on the matching progress is logged. If less than 1, then no progress information is logged. | `match`, `multimatch`, `filter` |
| `threads`      | `-1`       | The number of consumer threads *n* when processing data with respect to goals ``match``, ``filter`` and ``multimatch`` and also so during the update phase of the ``db`` goal. There is always one additional thread that reads and uncompresses a corresponding fastq or fasta file (so it is *n + 1* threads in total). When negative, the number of available processors *- 1* is used as *n*. When 0, then the corresponding goals run in single-threaded mode. | `db`, `match`, `multimatch`, `filter` |
| `countUniqueKMers`      | `true`       | If `true`, unique *k*-mers will be counted. This requires less than 5% of additional main memory.        |  `match`, `multimatch` |
| `writeDumpedFastq`   | `false`        | If `true`, then ``filter`` will also generate a fastq file `dumped_<fqfile>` with all reads not written to the corresponding filtered fastq file. |  `filter` |
| `writeFilteredFastq`   | `false`        | If `true`, then the goal `match` writes a filtered fastq file in the same way that the goal `filter` does. Moreover, Genestrip will write an output file `<fqfile>.out` in the [Kraken output format](https://ccb.jhu.edu/software/kraken/MANUAL.html#output-format) under `<base dir>/projects/<project_name>/krakenout` covering all filtered reads. | `match` |
| `matchWithKMerCounts`   | `false`        | Experimental: Counts how many times each unique *k*-mer has been detected.        |  `match`, `multimatch` |
| `maxKMerResCounts`   | `200`        | The number of the most frequent *k*-mers that will be reported, if `matchWithKMerCounts=true`.       |`match`, `multimatch` |
| `kMerSize`   | `31`        | The number of base pairs *k* for a *k*-mers. Changes to this values do *not* affect the memory usage of database. A value > 32 will cause collisions, i.e. leads to false positives for the `match` goal. | all |
| `useHttp` | `true` | Use http(s) to download data from NCBI. If ``false``, then Genestrip will try anonymous FTP instead (with login and password set to `anonymous`). | `db` |
| `refseqHttpBaseURL`   | `https://ftp.ncbi.nlm.nih.gov/refseq` | This [mirror](https://www.funet.fi/pub/mirrors/ftp.ncbi.nlm.nih.gov/refseq/) might be considered as an alternative. (No other mirror sites are known.) | `db` |
| `refseqFTPBaseURL`   | `ftp.ncbi.nih.gov`       |         | `db` |
| `taxHttpBaseURL`   | `https://ftp.ncbi.nlm.nih.gov`        | This base URL will be extended by the path `/pub/taxonomy/` in order to download the taxonomy file `taxdmp.zip`.        | all |
| `taxFTPBaseURL`   | `ftp.ncbi.nih.gov`        |         |  all |
| `ignoreMissingFastas`   | `true`        | If `true`, then a download of files from NCBI will not stop in case a file is missing on the [NCBI server](https://ftp.ncbi.nlm.nih.gov/).        |  `db` |
| `completeGenomesOnly`   | `false`        | If `true`, then only genomic accessions with the prefixes `AC`, `NC_`, `NZ_` will be considered when generating a database. Otherwise, all genomic accessions will be considered. See [RefSeq accession numbers and molecule types](https://www.ncbi.nlm.nih.gov/books/NBK21091/table/ch18.T.refseq_accession_numbers_and_mole/) for details.       |  `db` |
| `rankCompletionDepth` | *empty* | The rank up to which tax ids from `taxids.txt` will be completed by descendants of the taxonomy tree (the set rank included). If not set, the completion will traverse down to the lowest possible levels of the [taxonomy](https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip). Typical values could be `genus`, `species` or `strain`, but  all values used for assigning ranks in the taxonomy are possible. | `db` |
| `maxGenomesPerTaxid` | *unlimited* | The maximum number of genomes per tax id from the [RefSeq](https://ftp.ncbi.nlm.nih.gov/refseq/release/) to be included in the database. If negative, zero or not set, there is not limit. Note, that this is an important parameter to control database size, because in some cases, there are millions of genomic entries for a tax id such as for `573` (which does not even account for entries of its descendants). |  `db` |
| `useBloomFilterForMatch`   | `true`        | If `true` and if a bloom filter file for a database is present, it will be used during fastq file analysis (i.e. matching). Using the bloom filter tends to shorten matching time, if the most part of the reads cannot be classified because they contain *no* *k*-mers from the database. Otherwise, using the bloom filter might increase matching time by up to 30%. It also requires more main memory. | `match`, `multimatch` |
| `maxReadTaxErrorCount`   | `0.1`        | The absolute or relative maximum number of *k*-mers that do not have to be in the database for a read to be classified. If the number is above `maxReadTaxErrorCount`, then the read will not be classified. Otherwise the read will be classified in the same way as [done by Kraken](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2014-15-3-r46/figures/1).  If  `maxReadTaxErrorCount` is >= 1, then it is interpreted as an absolute number of *k*-mers. Otherwise (and so, if >= 0 and < 1), it is interpreted as the ratio between the *k*-mers not in the database and all *k*-mers of the read. | `match`, `multimatch` |
| `maxDust`   | `-1`        | When generating a database via the goal `db`, any low-complexity *k*-mer with too many repetitive sequences of base pairs may be omitted for storing. To do so, Genestrip employs a simple [genetic dust-filter](https://pubmed.ncbi.nlm.nih.gov/16796549/) for *k*-mers: It assigns a dust value *d* to each *k*-mer, and if *d* >  `maxDust`, then the *k*-mer will not be stored. Given a *k*-mer with *n* repeating base pairs of repeat length *k(1), ... k(n)* with *k(i) > 1*, then *d = fib(k(1)) + ... + fib(k(n))*, where *fib(k(i))* is the Fibonacci number of *k(i)*.  E.g., for the *8*-mer `TTTCGGTC`, we have *n = 2* with *k(1) = 3*, *k(2) = 2* and *d = fib(3) + fib(2) = 2 + 1 = 3*. For practical concerns `maxDust = 20` may be suitable. In this case, if *31*-mers were uniformly, randomly generated, then about 0.2 % of them would be omitted. If `maxDust = -1`, then dust-filtering is inactive.| `db` |
| `normalizedKMersFactor` | 1000000000 | A factor used to compute `normalized kmers` at read analysis time. | `match`, `multimatch` |
| `intialReadSizeBytes` | 4096 | The initial internal buffer size for reads in number of base pairs plus one. If longer reads occur, then internal buffer sizes will grow automatically. | `match`, `multimatch`, `filter` |
| `seqType` | `genomic` | Which type of sequence files to include from the RefSeq. Possible values are `genomic`, `rna` or `both`. RNA files from the RefSeq end with `rna.fna.gz`, whereas genomes end with `genomic.fna.gz`. | `db`, `index` |
| `classifyReads` | `true` | Whether to do read classification in the style of Kraken and KrakenUniq. Matching is faster without read classification and the columns `kmers`, `unique kmers` and `max contig length` in resulting CSV files are usually more conclusive anyways - in particular with respect to long reads. When read classification is off, the columns `reads` and `kmers from reads` will be 0 in resulting CSV files. | `match`, `multimatch` |
| `refSeqLimitForGenbankAccess` | `0` | Determines whether Genestrip should try to lookup genomic fasta files from Genbank, if the number of corresponding reference genomes from the RefSeq is below the given limit for a requested tax id. E.g. `refSeqLimitForGenbankAccess=1` would imply that Genbank is consulted if not a single reference genome is found in the RefSeq for a requested tax id. The default `refSeqLimitForGenbankAccess=0` essentially inactivates this feature. In addition, Genbank access is also influenced by the keys `fastaQualities` and `maxFromGenBank` (see below). | `db`, `index` |
| `fastaQualities` | `COMPLETE_LATEST` `,` `CHROMOSOME_LATEST` | Determines the allowed quality levels of fasta files from Genbank. The values must be comma-separated. The possible values are ordered from best to worst: `COMPLETE_LATEST`, `COMPLETE`, `CHROMOSOME_LATEST`, `CHROMOSOME`, `CONTIG_LATEST`, `CONTIG`, `LATEST`, `NONE`. If a corresponding value is included in the list, then a fasta file for a requested tax id on that quality level will be included, otherwise not (while also respecting the conditiions excerted via the keys `refSeqLimitForGenbankAccess` and `maxFromGenBank`). The quality levels are based on Genbank's [Assembly Summary File](https://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank.txt) (columns `version_status` and `assembly_level`). | `db`, `index` |
| `maxFromGenBank` | `1` | Determines the maximum number of fasta files used from Genbank per requested tax id. If the corresponding number of matching files exceeds `maxFromGenBank`, then then best ones according to `fastaQualities` will be retained to still match this maximum.  | `db`, `index` |


# Project properties

An optional configuration properties file `project.properties` or `Project.properties` may be put under the project folder `<base dir>/projects/<project_name>`.

The following entries are possible:
* `ignoreMissingFastas`,
* `completeGenomesOnly`,
* `rankCompletionDepth`,
* `maxGenomesPerTaxid`,
* `useBloomFilterForMatch`,
* `maxReadTaxErrorCount`,
* `maxDust`,
* `seqType` and
* `classifyReads`,
* `refSeqLimitForGenbankAccess`,
* `fastaQualities` and
* `maxFromGenBank`.

The use of the entries is the same as in the `config.properties` file. If given, an entry in `project.properties` overrides a corresponding entry from `config.properties` under this project.

# API-based usage

An API-based usage of the `match` goal is straight-forward: Please check out the test class `org.metagene.genestrip.APITest` in the folder `src/test/java` as a simple example.