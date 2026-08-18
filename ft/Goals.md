|Name|User Goal|Description|
|-|-|-|
|`ftgenall`|X|Generate the *k*-mer matching database and the filtering database according to Genestrip's `genall` and additionally the Genestrip-FT database with respect to the given project.|
|`ftsetup`||Create additional folders in `<base dir>/<project>` like `tex`.|
|`dendrogram`||Generate the dendrograms from *k*-mer intersection counts using agglomerative clustering.|
|`dendrolatex`|X|Generate LaTeX extracts for depicting the dendrograms from `dendrogram`.|
|`kmerindexbloom`||Store which *k*-mers belongs to which species for all *k*-mers under the rank genus (and potentially other ranks depending on configuration) in a Bloom filter.|
|`kmerindexsize`||Estimate how many entries the `kmerindexbloom` filter will hold, by reading the sequences once and sketching the (*k*-mer, genome) pairs with HyperLogLog. Only made when `kmerIndexSizing` asks for that estimate.|
|`intersectcount`||Count the number of joint *k*-mers between any two species per genus rank (and potentially other ranks depending on configuration).|
|`intersectcsv`|X|Save the number of joint *k*-mers between any two species per genus rank (and potentially other ranks depending on configuration) to CSV files along with resulting Jaccard-indices.|
|`branchhisto`||Compute, per tree node with refined *k*-mers, a histogram over the branching degrees of its *k*-mers, i.e. how many *k*-mers occur in exactly 1, 2, ... of the node's child subtrees (plus the trailing OTHER bucket as an additional child).|
|`branchhistocsv`||Write the branching-degree histograms from `branchhisto` to a CSV file, one row per tree node.|
|`branchhistorankcsv`||Aggregate the branching-degree histograms from `branchhisto` by taxonomic rank, one row per rank. Each node is turned into a relative branching-degree distribution (summing to one, excluding the OTHER column); for each branching degree from 1 up to 10, the mean, standard deviation, median and q1/q3 quartiles of that relative value across the rank's nodes are written as columns named `<degree>-avg`, `<degree>-stddev`, `<degree>-q1`, `<degree>-median`, `<degree>-q3`. Each row also starts with two summary blocks of the same five statistics: `childdeg-*` over the node's number of children in the database's taxonomy tree, and `kmerdeg-*` over the k-mer-weighted mean branching degree of the node (both excluding the OTHER column). Written to a CSV file.|
|`loadkmerindex`||Load the Bloom filter computed via the goal `kmerindexbloom`.|
|`storekmerindex`||Store the Bloom filter computed via the goal `kmerindexbloom`.|
|`ftupdatedb`||Update the database by integrating the refined taxonomy tree and reassigning *k*-mers under the genus ranks accordingly.|
|`ftdb`|X|Store the updated database.|
|`ftdbinfo`|X|Write information on the updated database content to a CSV file.|
|`loadftdb`||Load the updated database.|
|`allinonelatex`|X|Merge a project's LaTeX extracts from `dendrolatex` into one LaTeX document.|
|`ftmatchres`||Analyze fastq files according to Genestrip's `matchres` but with a Genestrip-FT database instead.|
|`ftmatch`|X|Analyze fastq files according to Genestrip's `match` but with a Genestrip-FT database instead.|
|`ftdb2fastq`|X|Generate fastq files according to Genestrip's `db2fastq` but from a Genestrip-FT database instead.|
|`ftclear`|X|Same as goal `clear`, but also clears `tex` the folder.|
|`ftsvgtaxtree`|X|Same as `svgtaxtree` but for an FT database.|
|`dbqualcounts`||Counts the *k*-mers per taxid directly from the underlying genomic files given a corresponding *k*-mer is in the database at all.|
|`ftqualcounts`||Same as `dbqualcounts` but for an FT database.|
|`dbquality`||Write the per-taxid quality metrics (tp, tp+fp, tp+fn, precision and recall) derived from `dbqualcounts` to a CSV file.|
|`ftquality`||Same as `dbquality` but for an FT database, i.e. derived from `ftqualcounts`.|
