|Name|User Goal|Description|
|-|-|-|
|`show`|X|Show user-related goals.|
|`showall`|X|Show user-related and most internal goals.|
|`genall`|X|Generate the *k*-mer matching database and also the filtering database with respect to the given project.|
|`clear`|X|Clear the folders `csv`, `db` and `krakenout`  of a project. This will delete all files the respective folders!|
|`db`|X|Generate the database for *k*-mer matching with respect to the given project.|
|`dbinfo`|X|Write information about a project's database content to a CSV file.|
|`tempdbinfo`||Write information about a project's temporary database content to a CSV file.|
|`db2fastq`|X|Generate fastq files from the database. A respective fastq file will contain all *k*-mers specifically associated with a single tax id from the database where each *k*-mer is represented by a read consisting of *k* bases. Respective fastq files will be stored in `<base dir>/projects/<project_name>/fastq` with the file name format `<project_name>_db2fastq_<taxid>.fastq.gz`. The command line option `tx` serves at selecting the corresponding tax ids for the fastq files to be generated. If the option is omitted, then fastq files for *all* tax ids from the database will be generated.|
|`index`|X|Generate a filtering database with respect to a given project.|
|`match`|X|Analyze fastq files as given by the `-f` or `-m` option. The resulting CSV file(s) will be stored in `<base dir>/projects/<project_name>/csv` unless specified otherwise via the `-r` option.|
|`matchlr`|X|Same as `match` but without doing read classification. This corresponds to the configuration setting `classifyReads=false`.|
|`filter`|X|Filter fastq files as given by the `-f` or `-m` option. The resulting filtered fastq file(s) `filtered_...` will be stored under `<base dir>/projects/<project_name>/fastq/` unless specified otherwise via the `-r` option.|
|`fasta2fastq`||Transform a fasta file or streaming resource to fastq file.|
|`matchres`||Analyze fastq files as given by the `-f` or `-m` option.|
|`matchreslr`||Same as `matchres` but without doing read classification.|
|`commonsetup`||Create data folders in `<base dir>/common` including `common` itself.|
|`setup`||Create data folders in `<base dir>/<project>`.|
|`taxdownload`||Download the taxonomy.|
|`taxtree`||Load the taxonomy into memory.|
|`taxnodes`||Compute the taxids for the project's database.|
|`refseqrelease`||Download the RefSeq release number.|
|`refseqcat`||Download the RefSeq release catalog files.|
|`checksummap`||Load MD5 checksums for RefSeq files into memory|
|`categories`||Load the RefSeq category names as requested for the project's database.|
|`refseqfna`||Download the genomic files from RefSeq for the requested categories.|
|`accmapsize`||Compute the number of RefSeq accession entries to be kept in memory.|
|`accmap`||Load the required RefSeq accession entries into memory.|
|`taxfromgenbank`||Determine for tax ids for which additional fasta files from Genbank should be downloaded.|
|`assemblydownload`||Download Genbank's assembly catalog file.|
|`fastasgenbank`||Determine the fasta files to be downloaded from Genbank.|
|`fastasgenbankdl`||Download the requested fasta files from Genbank.|
|`adddownloads`||Download additionally requested fasta files.|
|`addfastas`||Associate additional fasta files with tax ids in memory.|
|`fillsize`||Precompute the number of *k*-mers for the project's database.|
|`tempindex`||Fill the temporary bloom index with *k*-mers.|
|`filldb`||Fill the database with *k*-mers.|
|`tempdb`||Save temporary database.|
|`loadtempdb`||Load the temporary database into memory.|
|`updatedb`||Update the temporary database with regard to least commonn ancestors of *k*-mers.|
|`loaddb`||Load the matching database into memory.|
|`fillindex`||Fill the filtering database with *k*-mers.|
|`loadindex`||Load the filtering database into memory.|
|`db2fastqtaxids`||Determine tax ids for the goal `db2fast` from the command line argument.|
|`fastqmap`||Determine the fastq mapping from command line arguments.|
|`fastqmaptransform`||Transform URLs of fastq files to be downloaded to local paths.|
|`fastamap`||Determine the fasta mapping from command line arguments.|
|`fastamaptransform`||Transform URLs of fasta files to be downloaded to local paths.|
|`fastqdownload`||Download fastq files given via URLs as requested.|
|`fastadownload`||Download fasta files given via URLs as requested.|
|`krakencount`||For internal use (to invoke kraken and count results).|
|`krakenres`||For internal use (to write kraken results to a file).|
|`dbdownload`||Download and install a project's database via a given URL.|
