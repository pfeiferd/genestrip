
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

**Genestrip**: Efficient read classification, filtering and *k*-mer counting for selected groups of species 
===============================================

## Introduction

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

## License

[Genestrip is free for non-commercial use.](./LICENSE.txt) Please contact [daniel.pfeifer@progotec.de](mailto:daniel.pfeifer@progotec.de) if you are interested in a commercial license.

## Building and installing

Genestrip is structured as a standard [Maven 2 or 3](https://maven.apache.org/) project and is compatible with the [JDK](https://jdk.java.net/) 11 or higher.[^1]

To build it, `cd` to the installation directory `genestrip`. Given a matching Maven and JDK installation, `mvn install` will compile and install the Genestrip program library. Afterwards a self-contained and executable `genestrip.jar` file will be stored under `./lib`. 

Since version 0.5, Genestrip is also available on [Maven Central](https://repo1.maven.org/maven2/org/genestrip/genestrip/).
Here  is the dependency:
```
<dependency>
	<groupId>org.genestrip</groupId>
	<artifactId>genestrip</artifactId>
	<version>2.1</version>
</dependency>
```
You may check for higher versions and update the dependency accordingly...

### Self-contained Windows executable

The Maven command `mvn -P winexe package` creates a self-contained Windows executable `bin/genestrip.exe` along with a reduced JRE under `lib/jre`. 
You need a [JDK](https://jdk.java.net/) 11 installation or higher on Windows for this command to succeed.
Moreover it creates a Zip folder `target/genestrip-${version}-windows-x64.zip` that contains the same files with less than 40 MB of disk space.

When extracted Genestrip can be executed via `bin\genestrip.exe` (on a Windows x86, 64 bit architecture). It will search for the above JRE for 
execution under the relative location `..\lib\jre`. So the folders `bin` and `lib` should be kept together accordingly. **These few files is all it takes to run Genestrip under Windows!** There is no additional JRE necessary.

For convenience, the ready-made Zip folder `genestrip-${version}-windows-x64.zip` is also [publicly available for download on genestrip.it.hs-heilbronn.de](https://genestrip.it.hs-heilbronn.de/files/bin).

### Runnning the JUnit tests

The Maven command `mvn -P prerelease install` runs all the JUnit tests for Genestrip and more. It will take time...

[^1]: Counter to common belief, Java can well be used for such high performance applications when using its programming facilities the right way.

## Genestrip databases

### Generating the sample database

The Genestrip installation holds additional folders that include the sample project `human_virus`. After building Genestrip, you may call
`sh ./bin/genestrip.sh human_virus dbinfo`
in order to generate the `human_virus` database and create a CSV file with basic information about the database content.

Genestrip follows a goal-oriented approach in order to create any result files (in similarity to [make](https://pubs.opengroup.org/onlinepubs/9699919799/utilities/make.html)). So, when generating the `human_virus` database, Genestrip will
1. download the [taxonomy](https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip) and unzip it to `./data/common`,
1. download the [RefSeq release catalog](https://ftp.ncbi.nlm.nih.gov/refseq/release/release-catalog/) to `./data/common/refseq`,
1. download [all virus related RefSeq fna files](https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/) to `./data/commmon/refseq` (which is currently just a single file),
1. perform several follow-up goals, until the database file `human_virus_db.zip` is finally made under `./data/projects/human_virus/db`, 
1. create a CSV file `human_virus_dbinfo.csv` under `./data/projects/human_virus/csv`, which contains basic information about the database, i.e. the number of *k*-mers stored per tax id.

The generated database comprises *k*-mers for all viruses according to the tax id file `./data/project/human_virus/taxids.txt`.


### Generating your own database

Generating your own database is straight-forward:
1. Create a project folder `<project_name>` under `./data/projects/`. This is the place for all subsequent files that specify the content of a database to be generated. It is also the core name of a related database.
1. Create a text file `./data/projects/<project_name>/taxids.txt` with one tax id per line. The database will *only* contain *k*-mers from genomes that belong to tax ids referenced in the file `taxids.txt`. If `taxids.txt` contains a tax id from a non-leaf node of the [taxonomy](https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip), then all subordinate tax ids and *k*-mers from respective [RefSeq](https://ftp.ncbi.nlm.nih.gov/refseq/release/) genomes will be included in the database as well.
1. Create a text file `./data/projects/<project_name>/categories.txt`. This file tells Genestrip, which categories of organisms should be *considered* when generating the database and when determining a least common ancestor tax id for *k*-mers in the finalized database. You also have to ensure that the categories from `categories.txt` comprise all of the tax ids from `taxids.txt`: Only genomes of comprised tax ids will be used to generate the database.
1. The following categories are allowed, and originate from the [RefSeq](https://ftp.ncbi.nlm.nih.gov/refseq/release/):
`archaea`, `bacteria`, `complete`, `fungi`, `invertebrate`, `mitochondrion`, `other`, `plant`, `plasmid`, `plastid`, `protozoa`, `vertebrate_mammalian`, `vertebrate_other` and `viral`. You may enter one category per line.
The *entire* set of genomes that belong to a category referenced in `categories.txt` will be used to determine the least common ancestors tax ids of the *k*-mers stored in the finalized database. The more categories you reference in `categories.txt`, the more genomes files will have to be downloaded from the [RefSeq](https://ftp.ncbi.nlm.nih.gov/refseq/release/) and the longer the database generation process will take. E.g., if only the category `viral` is included, it should just take a few minutes, but with `bacteria` included, it may typically take about a day or more.
So you should choose a suitable set of categories to balance completeness of considered genomes with efficiency of the database generation process.

1. Start the database generation process via the command `sh ./bin/genestrip.sh <project_name> dbinfo`. 
This will also create a CSV file `./data/projects/<project_name>/csv/<project_name>_dbinfo.csv` with 
basic information about the database, i.e. the number of *k*-mers stored per tax id. 
The path of the database will be `./data/projects/<project_name>/db/<project_name>_db.zip`. 
The database file `<project_name>_db.zip` is self-contained and includes all necessary taxonomic information for its 
later use via the goals `match` and `matchlr`.

In general, Genestrip organizes a project folder `./data/projects/<project_name>` by means of the following sub-folders:
* `csv` is where analysis results of fastq files will be stored (by default).
* `db` is where Genestrip puts the generated database. If your focus is on filtering fastq files, Genestrip can create a specialized,  filtering database named `<project_name>_index.ser.gz` that will be put there too. Moreover, the intermediate database `<project_name>_tempdb.zip` will be put there temporarily as part of the database generation process. 
* `fasta` may be used to store *additional* genome files (not part of in the [RefSeq](https://ftp.ncbi.nlm.nih.gov/refseq/release/)) that should be considered for the database generation process (see Section [Additional fasta files](#additional-fasta-files)).
* `fastq` is where Genestrip will store filtered (or generated) fastq files. You may also put you own fastq files to be analyzed there (but they can be read from any other path too).
* `genbank` is where genomic files from Genbank will be downloaded to if needed during the project's database generation process.
* `krakenout` is for output files in the style of [Kraken](https://ccb.jhu.edu/software/kraken/MANUAL.html#output-format). They may optionally be generated when analyzing fastq files.
* `log` is reserved for future use in conjunction with a server-side integration of Genestrip.

A database covering more species may require more memory - especially while generating the database. This can be addressed by adjusting the script `genestrip.sh` where the argument `-Xmx32g` sets the maximum heap space of the Java Virtual Machine to 32GB. E.g. to double it, simple replace `32g` by `64g`.

### Some preconfigured and ready-made databases

There is a separate project [Genestrip-DB](https://github.com/pfeiferd/genestrip-db) on GitHub that covers [several databases](https://github.com/pfeiferd/genestrip-db/blob/main/README.md#the-databases) and offers a [corresponding download from genestrip.it.hs-heilbronn.de](https://genestrip.it.hs-heilbronn.de/files/data). 

## Analysis with Genestrip

### Analyzing fastq files by matching *k*-mers of reads

Genestrip's main purpose is to analyze reads from fastq files and count the contained *k*-mers per tax id according to a previously generated database. As an example, Genestrip comes with a small fastq-file `sample.fastq.gz` in `./data/projects/human_virus/fastq`. To start the matching process for it, run
```
sh ./bin/genestrip.sh human_virus match -f ./data/projects/human_virus/fastq/sample.fastq.gz
```
The resulting CSV file will be named `human_virus_match_sample.csv` under `./data/projects/human_virus/csv`. The same applies to your own projects under `./data/projects`.

These are a few lines of its contents along with the header line:
```
pos;name;rank;taxid;reads;kmers from reads;kmers;unique kmers;contigs;average contig length;max contig length;reads >=1 kmer;read bps;avg. read length;db coverage;exp. unique kmers;unique kmers / exp.;db kmers;parent taxid;norm. reads;norm. kmers;norm. reads bps;norm. read >=1 kmer;norm. reads kmers;acc. reads;acc. norm. reads;acc. kmers;acc. norm. kmers;acc. reads bps;acc. norm. reads bps;acc. read >=1 kmer;acc. norm. read >=1 kmer;acc. reads kmers;acc. norm. reads kmers;max contig desc.;
0;TOTAL;;;6565;0;461305;0;0;;0;0;0;;;;;3951346;;;;;;;;;;;;;;;;;;
1;root;no rank;1;0;0;0;0;0;;0;0;0;;;0.0;;0;;;;;;;6565;23.985741934304116;50829;568.1616183493833;658255;2401.884786157381;10050;170.9186006856684;50814;236.08784005861423;;
2;Viruses;superkingdom;10239;1155;11850;28388;73;3616;7.850663716814159;68;2988;115668;100.14545454545454;0.9125;80.0;0.9125;80;1;14.4375;354.85;1445.85;37.35;148.125;6565;23.985741934304116;50829;568.1616183493833;658255;2401.884786157381;10050;170.9186006856684;50814;236.08784005861423;@NS500362:54:HT523BGX2:4:13606:8875:6204 2:N:0:2;
3;Riboviria;clade;2559587;0;0;0;0;0;;0;0;0;;;0.0;;0;10239;;;;;;3442;9.473678636821118;16166;202.3210751319162;345107;948.5733099243974;4822;123.68058656389968;22441;87.4796036347947;;
4;Orthornavirae;kingdom;2732396;25;277;240;4;237;1.0126582278481013;2;236;2502;100.08;1.0;4.0;1.0;4;2559587;6.25;60.0;625.5;59.0;69.25;3442;9.473678636821118;16166;202.3210751319162;345107;948.5733099243974;4822;123.68058656389968;22441;87.4796036347947;@NS500362:54:HT523BGX2:4:11405:19113:11578 2:N:0:2;
5;Negarnaviricota;phylum;2497569;0;0;0;0;0;;0;0;0;;;0.0;;0;2732396;;;;;;50;0.003993291270665282;108;0.00862550914463701;5018;0.40076671192396773;50;0.003993291270665282;108;0.00862550914463701;;
6;Polyploviricotina;subphylum;2497571;0;0;0;0;0;;0;0;0;;;0.0;;0;2497569;;;;;;50;0.003993291270665282;108;0.00862550914463701;5018;0.40076671192396773;50;0.003993291270665282;108;0.00862550914463701;;
7;Ellioviricetes;class;2497576;0;0;0;0;0;;0;0;0;;;0.0;;0;2497571;;;;;;50;0.003993291270665282;108;0.00862550914463701;5018;0.40076671192396773;50;0.003993291270665282;108;0.00862550914463701;;
8;Bunyavirales;order;1980410;0;0;0;0;0;;0;0;0;;;0.0;;0;2497576;;;;;;50;0.003993291270665282;108;0.00862550914463701;5018;0.40076671192396773;50;0.003993291270665282;108;0.00862550914463701;;
9;Phenuiviridae;family;1980418;0;0;0;0;0;;0;0;0;;;0.0;;0;1980410;;;;;;50;0.003993291270665282;108;0.00862550914463701;5018;0.40076671192396773;50;0.003993291270665282;108;0.00862550914463701;;
10;Phlebovirus;genus;11584;0;0;0;0;0;;0;0;0;;0.0;0.0;;41;1980418;0.0;0.0;0.0;0.0;0.0;50;0.003993291270665282;108;0.00862550914463701;5018;0.40076671192396773;50;0.003993291270665282;108;0.00862550914463701;;
11;Phlebovirus toroense;species;3052684;0;0;0;0;0;;0;0;0;;;0.0;;0;11584;;;;;;50;0.003993291270665282;108;0.00862550914463701;5018;0.40076671192396773;50;0.003993291270665282;108;0.00862550914463701;;
12;Punta Toro virus;no rank;11587;50;108;108;15;52;2.076923076923077;9;50;5018;100.36;0.0011979873811995847;107.53983475448534;0.13948319740536305;12521;3052684;0.003993291270665282;0.00862550914463701;0.40076671192396773;0.003993291270665282;0.00862550914463701;50;0.003993291270665282;108;0.00862550914463701;5018;0.40076671192396773;50;0.003993291270665282;108;0.00862550914463701;@NS500362:54:HT523BGX2:4:21511:11527:10410 2:N:0:2;
```
[**This is a list of all CSV columns.**](CSVColumns.md)

The frequencies from the column `max kmer counts` can be used to build frequency graph's for *k*-mers as shown below. The frequency graphs help to further assess the validity of analysis results.

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

### Reliability of results

We cannot guarantee for any results returned by Genestrip. Use this software at your own risk. **Important: It is by no means meant to be used for any medical purposes** and it is purely experimental in nature.

Despite of the these limitations, we tested the Genestrip in the following ways:
* Critical code is covered by functional tests (using [JUnit](https://junit.org)).
* Results have been evaluated using various benchmark datasets.
* We applied Genestrip to several real-world fastq files, 
where the findings matched the expectations. For example, we correctly recovered bacteria in ticks from fastq files as described in [this publication](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10328957/).

### Filtering fastq files

Genestrip can be used to *filter* fastq files via *k*-mers from a previously generated database. As an example, one may also use the fastq file `sample.fastq.gz` from `./data/projects/human_virus/fastq`. To start the filtering process, run
````
sh ./bin/genestrip.sh human_virus filter -f ./data/projects/human_virus/fastq/sample.fastq.gz
````
First, the command creates a filtering database file `human_virus_index.ser.gz`, if not yet present, under `./data/projects/<project_name>/db`. 

The resulting filtered fastq file named `human_virus_filtered_sample.fastq.gz` will be put under `./data/projects/human_virus/fastq`. It holds just the reads which contain at least one *k*-mer from the `human_virus` database. Note that, in case of the specific sample fastq file, the filtered fastq file has about the same size as the original file because the sample had *originally been created to just contain human virus-related reads*.

A filtering database is typically smaller than a database required for *k*-mer matching but the former can only be used for filtering. So, if you are tight on resources and your focus is on filtering, you may prefer using a filtering database. Also, the filtering process is  faster than the *k*-mer matching process.

## Technical documentation

### Usage and goals

The usage of Genestrip:
```
usage: genestrip [options] <project> [<goal1> <goal2>...]
 -C <key>=<value>           To set Genestrip configuration paramaters via
                            the command line.
 -d <base dir>              Base directory for all data files. The default
                            is './data'.
 -db <database>             Path to filtering or matching database for the
                            goals 'filter' or 'match', 'matchlr', 'dbinfo'
                            and 'db2fastq' for use without project
                            context.
 -f <fqfile1,fqfile2,...>   Input fastq files as paths or URLs to be
                            processed via the goals 'filter', 'match' or
                            'matchlr'. When a URL is given, the fastq file
                            will not be downloaded but data streaming will
                            be applied unless '-l' or '-ll' is given.
 -k <key>                   Key used as a prefix for naming result files
                            in conjuntion with '-f'.
 -l                         Download fastqs from URLs to '<base
                            dir>/projects/<project name>/fastq' instead of
                            streaming them for the goals 'filter', 'match'
                            and 'matchlr'.
 -ll                        Download fastqs from URLs to '<base
                            dir>/fastq' instead of streaming them for the
                            goals 'filter', 'match' and 'matchlr'.
 -m <fqmap>                 Mapping file with a list of fastq files to be
                            processed via the goals 'filter', 'match' or
                            'matchlr'. Each line of the file must have the
                            format '<key> <URL or path to fastq file>'.
 -r <path>                  Common store folder for filtered fastq files
                            and result files created via the goals
                            'filter', 'match' or 'matchlr'. The defaults
                            are '<base dir>/projects/<project name>/fastq'
                            and '<base dir>/projects/<project name>/csv',
                            respectively.
 -t <target>                Generation target ('make', 'clean' or
                            'cleanall'). The default is 'make'.
 -tx <taxids>               List of tax ids separated by ',' (but no
                            blanks) for the goal 'db2fastq'. A tax id may
                            have the suffix '+', which means that
                            taxonomic descendants from the project's
                            database will be included.
 -v                         Print version.                            
```

Genestrip follows a goal-oriented approach in order to create any result file (in similarity to  [make](https://pubs.opengroup.org/onlinepubs/9699919799/utilities/make.html)). Goals are executed in a lazy manner, i.e. a file is only (re-)generated, if it is missing at its designated place in the `<base dir>` folder or any of its subfolders.

Most user-related goals are project-oriented. This means that a database project must exist for the goal to be executable. The goals `match`, `matchlr` and `filter` refer to fastq files. This means that they can only be executed if fastq files are given via the `-f` or `-m` option. 

Genestrip supports matching of *k*-mers from multiple fastq files in batches: For this purpose you may create a mapping file with one line per fastq file to be processed. Each line should have the following form:
```
<key> <path_or_URL_to_fastq_file>
```
If several fastq files in the mapping file have the same `<key>`, then the matching results of these files will be merged. 
A resulting CSV file named `<project_name>_imatch_<key>.csv` will be put under `<base dir>/projects/<project_name>/csv` unless specified otherwise via the `-r` option. 
If `<path_or_URL_to_fastq_file>` is a file name without a path prefix, the file is assumed to be located in `<base dir>/projects/<project_name>/fastq` or otherwise in `<base dir>/fastq`.

[Glob patterns](https://docs.oracle.com/javase/tutorial/essential/io/fileOps.html#glob) may be used for all file names (i.e. the last component of a file path) in `<path_or_URL_to_fastq_file>`:
If the glob pattern matches a file from one of the above directories, then *all* matching files from that directory will be added to the joint analysis. 

If `<path_or_URL_to_fastq_file>` is a URL then the fastq file will be streamed or downloaded depending on the command line options `-l` or `-ll` (see Section [Reading, streaming and downloading fastq files](#reading-streaming-and-downloading-fastq-files)).

Fastq files and fasta files may be g-zipped or not. Genestrip will automatically recognize g-zipped files via the suffixes `.gz` and `.gzip`.

Some named goals are for internal purposes only. In principle, they could be run directly by users but rather serve the generation process of databases or they exist for experimental reasons.

[**This is a list of all goals.**](Goals.md)

Many goals depend on other goals. E.g., the `dbinfo` goal requires the corresponding database to exist and so, it will trigger the execution of the ``db`` goal in case the corresponding database is missing and so on. 

The following goal graph depicts the goals' dependencies (without the trivial goal `setup` because of too many edges).
`o:...` are goals whose result is an object in memory. `f:...` are file goals that produce one or more files as a result. `d:...` are download goals - they download files from a source and store them locally for further processing.
<p align="center">
  <img src="GoalGraph.svg" width="1400"/>
</p>

### Reading, streaming and downloading fastq files

Regarding analysis, fastq files will processed in one of the following ways:
1. Reading from the file system: This happens if a *file path* is given after the `-f` option.
1. *Streaming* from the network: This happens if a *URL* is given after the `-f` option. (The corresponding fastq file will not be downloaded to the file system then.)
1. *Downloading* from the network: This happens if a *URL* is given after the `-f` option *and if the option `-l` or `-ll`* is given, respectively. The corresponding fastq file be downloaded if not yet present. Afterwards, the downloaded file will be processed from the file system.

A comma-separated list of file paths or URLs or both may be put after `-f` without blanks. In this case, related fastq files will be analyzed together and the results will be merged.

Unless a key is given via `-k`, a resulting CSV-file will be named after the first file path or URL as given via `-f`.
Otherwise the name of key will be used for it. E.g.:
```
sh ./bin/genestrip.sh human_virus match -k mykey -f ./data/projects/human_virus/fastq/sample.fastq.gz
```
will result in the CSV file `./data/projects/human_virus/csv/human_virus_match_mykey.csv`.

In order to run the goals `match` or `matchlr` a project context is not strictly needed, when a database file is given via the option `-db`. 
E.g., the following command works with an arbitrary project name (here `someprojectname`):
```
sh ./bin/genestrip.sh someprojectname match -r . -k mysample -f ./data/projects/human_virus/fastq/sample.fastq.gz -db ./data/projects/human_virus/db/human_virus_db.zip
```
and produces the result file `./someprojectname_match_mysample.csv`.

### Targets

Genestrip supports three targets for each goal, namely `make`, `clean` and `cleanall`.

- `make` is the default target. It executes a given goal as explained before and creates the files associated with the goal in a lazy manner. If the corresponding files already exist, then `make` does nothing.
- `clean` *deletes* all files associated with the given goal but it does not delete any files of goals that the given goal depends on.
- `cleanall` does same as `clean`, but it *also recursively deletes* any files of goals that the given goal depends on.

### Configuration parameters

An optional configuration properties file `config.properties` may be put under `<base dir>`.
Entries per line should have the form
```
<key>=<value>
```
[**This is a list of all configuration parameters**](ConfigParams.md)

An optional configuration properties file `config.properties` may also be put under the project folder `<base dir>/projects/<project_name>`. Configuration entries from the project level override entries from the level `<base dir>`.

Moreover, configuration parameters may be set on the command line like this:
```
-C <key>=<value>
```
They have the highest priority.

### Additional fasta files

The following configuration files exclusively affect the generation of databases via the goals `db` and `index`.

#### Manually adding fasta files

In some cases you may want to add *k*-mers of genomes to your database, where the genomes are not part of the [RefSeq](https://ftp.ncbi.nlm.nih.gov/refseq/release/) (or Genbank). Genestrip supports this via an optional text file `additional.txt` under `<base dir>/projects/<project_name>`.
The file should contain one line for each additional genome file. (The genome file must be in fasta format and may be g-zipped or not.)
The line format is
```
<taxid> <path_to_fasta_file>
```
where `<taxid>` is the (unique) tax id associated with the file's genomic data. (Multiple tax ids per fasta file are not supported in this context.) 
If `<path_to_fastq_file>` is a file name without a path prefix, then the file is assumed to be located in `<base dir>/projects/<project_name>/fasta`. If not found there,
the directory `<base dir>/common/fasta` will be checked as a secondary location.

This adding of fasta files can also be used to *just* correct the least common ancestor of *k*-mers in the resulting database since the added fasta files will be automatically used during the update phase of the ``db`` goal. E.g., to correct the least common ancestor of *k*-mers occurring in a purely `protozoa`n database *but also* in the human genome, one may simple add
```
9606 <path_to_human_genome_fasta_file>
```
to `additional.txt` where `<path_to_human_genome_fasta_file>` points to the [human genome fasta file](https://www.ncbi.nlm.nih.gov/datasets/taxonomy/9606/). (Note that for this purpose, ``9606`` *must not* occur in `taxids.txt`, since otherwise all *k*-mers from the human genome would be included in a `protozoa`n database.)

#### Automated download of additional fasta files per project

The manually adding of fasta files as described above involves the manual download and positioning in the local file system.
To automate the download, *an extended format for entries in `additional.txt` is also possible*:
```
<taxid> <path_to_fasta_file> <URL_to_fasta_file> [<md5 fingerprint>]
```
The `<md5 fingerprint>` is optional and if given, will be used to ensure consistency of the downloaded file.
E.g., to automate the download of the human genome, the following entry will suffice:
```
9606 GCA_000001405.29_GRCh38.p14_genomic.fna.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.29_GRCh38.p14/GCA_000001405.29_GRCh38.p14_genomic.fna.gz 22907ad69ddfae66071bf9cf99b3e8de
```
The corresponding download file will be stored under `<base dir>/projects/<project_name>/fasta`.

#### Automated download of additional fasta files across projects

The automated download from the previous section is unsuitable, if the fasta file is large and is needed in several projects, as it will be downloaded and stored once per project.
To enable an automated download *across projects*, the file `<base dir>/common/fasta/downloads.txt` may be created with the line format:
```
<fasta_file_name> <URL_to_fasta_file> [<md5 fingerprint>]
```
So in this case, if a file named `<fasta_file_name>` cannot be found in `<base dir>/projects/<project_name>/fasta` or `<base dir>/common/fasta`, it will be downloaded
and stored as `<base dir>/common/fasta/<fasta_file_name>`. Afterwards it will be available with regard to entries like
```
<taxid> <fasta_file_name>
```
from a project's `additional.txt` file.

### API-based usage

An API-based invocation of the goals `match` and `filter` is straight-forward: Please check out the test class [`org.metagene.genestrip.APITest`](./src/test/java/org/metagene/genestrip/APITest.java) in the folder `src/test/java` as a code example.