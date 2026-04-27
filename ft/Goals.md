|Name|User Goal|Description|
|-|-|-|
|`ftsetup`||Create additional folders in `<base dir>/<project>` like `tex`.|
|`dendrogram`||Generate the dendrograms from *k*-mer intersection counts using agglomerative clustering.|
|`dendrolatex`|X|Generate LaTeX extracts for depicting the dendrograms from `dendrogram`.|
|`kmerindexbloom`||Store which *k*-mers belongs to which species for all *k*-mers under the rank genus (and potentially other ranks depending on configuration) in a Bloom filter.|
|`intersectcount`||Count the number of joint *k*-mers between any two species per genus rank (and potentially other ranks depending on configuration).|
|`intersectcsv`|X|Save the number of joint *k*-mers between any two species per genus rank (and potentially other ranks depending on configuration) to CSV files along with resulting Jaccard-indices.|
|`loadkmerindex`||Load the Bloom filter computed via the goal `kmerindexbloom`.|
|`storekmerindex`||Store the Bloom filter computed via the goal `kmerindexbloom`.|
|`updatestore`||Update the database by integrating the refined taxonomy tree and reassigning *k*-mers under the genus ranks accordingly.|
|`ftdb`|X|Store the updated database.|
|`ftdbinfo`|X|Write information on the updated database content to a CSV file.|
|`loadftdb`||Load the updated database.|
|`allinonelatex`|X|Merge a project's LaTeX extracts from `dendrolatex` into one LaTeX document.|
|`ftmatchres`||Analyze fastq files according to Genestrip's `matchres` but with a Genestrip-FT database instead.|
|`ftmatch`|X|Analyze fastq files according to Genestrip's `match` but with a Genestrip-FT database instead.|
|`ftdb2fastq`|X|Generate fastq files according to Genestrip's `db2fastq` but from a Genestrip-FT database instead.|
|`ftclear`|X|Same as goal `clear`, but also clears `tex` the folder.|
|`ftsvgtaxtree`|X|Same as `svgtaxtree` but for an FT database.|
|`tax2kmers`||Counts the *k*-mers per taxid directly from the underlying genomic files given a corresponding *k*-mer is in the database at all.|
|`fttax2kmers`||Same as `tax2kmers` but for FT database.|
|`ftcons`|X|TODO|
|`dbcons`|X|TODO|
