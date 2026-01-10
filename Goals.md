|Name|User Goal|Rec. Clean|Description|
|-|-|-|-|
|`show`|X|X|Show user-related goals.|
|`showall`|X|X|Show user-related and most internal goals.|
|`genall`|X|X|Generate the *k*-mer matching database and also the filtering database with respect to the given project.|
|`clear`|X|X|Clear the folders `csv`, `db` and `krakenout`  of a project. This will delete all files the respective folders!|
|`db`|X|X|Generate the database for *k*-mer matching with respect to the given project.|
|`dbinfo`|X|X|Write information on a project's database content to a CSV file.|
|`tempdbinfo`||X|Write information about a project's temporary database content to a CSV file.|
|`db2fastq`|X|X|Generate fastq files from the database. A respective fastq file will contain all *k*-mers specifically associated with a single tax id from the database where each *k*-mer is represented by a read consisting of *k* bases. Respective fastq files will be stored in `<base dir>/projects/<project_name>/fastq` with the file name format `<project_name>_db2fastq_<taxid>.fastq.gz`. The command line option `tx` serves at selecting the corresponding tax ids for the fastq files to be generated. If the option is omitted, then fastq files for *all* tax ids from the database will be generated.|
|`index`|X|X|Generate a filtering database with respect to a given project.|
|`match`|X|X|Analyze fastq files as given by the `-f` or `-m` option. The resulting CSV file(s) will be stored in `<base dir>/projects/<project_name>/csv` unless specified otherwise via the `-r` option.|
|`matchlr`|X|X|Same as `match` but without doing read classification. This corresponds to the configuration setting `classifyReads=false`.|
|`filter`|X|X|Filter fastq files as given by the `-f` or `-m` option. The resulting filtered fastq file(s) `filtered_...` will be stored under `<base dir>/projects/<project_name>/fastq/` unless specified otherwise via the `-r` option.|
|`extract`|X|X|Extract reads from fastq files based on matching descriptors. See also config key `extractKey`.|
|`fasta2fastq`||X|Transform a fasta file or streaming resource to fastq file.|
|`matchres`||X|Analyze fastq files as given by the `-f` or `-m` option.|
|`matchreslr`||X|Same as `matchres` but without doing read classification.|
|`commonsetup`||X|Create data folders in `<base dir>/common` including `common` itself.|
|`setup`||X|Create data folders in `<base dir>/<project>`.|
|`taxdownload`|||Download the taxonomy.|
|`taxtree`||X|Load the taxonomy into memory.|
|`taxnodes`||X|Compute the taxids for the project's database.|
|`refseqrelease`||X|Download the RefSeq release number.|
|`refseqcat`|||Download the RefSeq release catalog files.|
|`checksummap`|||Load MD5 checksums for RefSeq files into memory|
|`categories`||X|Load the RefSeq category names as requested for the project's database.|
|`refseqfna`|||Download the genomic files from RefSeq for the requested categories.|
|`accmapsize`||X|Compute the number of RefSeq accession entries to be kept in memory.|
|`accmap`||X|Load the required RefSeq accession entries into memory.|
|`taxfromgenbank`||X|Determine for tax ids for which additional fasta files from Genbank should be downloaded.|
|`assemblydownload`|||Download Genbank's assembly catalog file.|
|`fastasgenbank`||X|Determine the fasta files to be downloaded from Genbank.|
|`fastasgenbankdl`|||Download the requested fasta files from Genbank.|
|`adddownloads`|||Download additionally requested fasta files.|
|`addfastas`||X|Associate additional fasta files with tax ids in memory.|
|`fillsize`||X|Precompute the number of *k*-mers for the project's database.|
|`tempindex`||X|Fill the temporary bloom index with *k*-mers.|
|`filldb`||X|Fill the database with *k*-mers.|
|`tempdb`||X|Save temporary database.|
|`loadtempdb`||X|Load the temporary database into memory.|
|`updatedb`||X|Update the temporary database with regard to least commonn ancestors of *k*-mers.|
|`loaddb`||X|Load the matching database into memory.|
|`fillindex`||X|Fill the filtering database with *k*-mers.|
|`loadindex`||X|Load the filtering database into memory.|
|`db2fastqtaxids`||X|Determine tax ids for the goal `db2fast` from the command line argument.|
|`fastqmap`||X|Determine the fastq mapping from command line arguments.|
|`fastqmaptransform`||X|Transform URLs of fastq files to be downloaded to local paths.|
|`fastamap`||X|Determine the fasta mapping from command line arguments.|
|`fastamaptransform`||X|Transform URLs of fasta files to be downloaded to local paths.|
|`fastqdownload`|||Download fastq files given via URLs as requested.|
|`fastadownload`|||Download fasta files given via URLs as requested.|
|`krakencount`||X|For internal use (to invoke kraken and count results).|
|`krakenres`||X|For internal use (to write kraken results to a file).|
|`dbdownload`|||Download and install a project's database via a given URL.|
|`checkrefseqrnum`||X|Check whether the downloaded RefSeq release is equal to the current release on the download server.|
