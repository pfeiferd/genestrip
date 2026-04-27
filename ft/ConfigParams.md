|Name|Type|Value Range|Default|Description|For Goals|
|-|-|-|-|-|-|
|`clusterMethod`|nominal|`SINGLE_LINKAGE`, `COMPLETE_LINKAGE`, `UPGMA`, `WPGMA`|`SINGLE_LINKAGE`|The cluster distance to be used when performing agglomerative clustering.|`dendrogram`|
|`withDescendantCounts`|boolean||`false`|Whether to include the *k*-mer counts of *all* descendents for any two considered species in the denominator of the Jaccard-index. If not, only the *k*-mer counts right for the two considered species are used.|`intersectcount`|
|`turnLatex`|boolean||`true`|Whether the dendrogram in LaTeX has the species names aligned horizontally (with the entired diagram turned) or not.|`dendrolatex`|
|`xFactorLatex`|double|[0.0, 1.7976931348623157E308]|`1.0`|The factory for stretching the dendrogram in TikZ's native *x* coordinate.|`dendrolatex`|
|`yFactorLatex`|double|[0.0, 1.7976931348623157E308]|`8.0`|The factory for stretching the dendrogram in TikZ's native *y* coordinate.|`dendrolatex`|
|`tikzScaleFactor`|double|[0.0, 1.7976931348623157E308]|`1.0`|The factor for `scale` in the 'tikzpicture' environment of a LaTex dendrogram.|`dendrolatex`|
|`simLogScaling`|boolean||`false`|Whether to do logarithmic scaling of the similarity `sim` in dendrograms (via `1 - log (sim) / log (min_sim)`.|`dendrolatex`|
|`ftBloomFilterFpp`|double|[0.0, 1.0]|`0.001`|False positive probability of the Bloom filter used in Genestrip-FT.|`kmerindexbloom`|
|`allInOneChunkSize`|int|[1, 2147483647]|`50`|Maximum number of dendrograms put in one LaTeX file via the goal `allinonelatex`.|all|
|`refinementPositions`|comma-separated list of values of `<rank>`, `<taxid>` or else `*` which means all ranks and taxids are included|<rank> as subset of `cellular root`, `acellular root`, `superkingdom`, `domain`, `realm`, `kingdom`, `phylum`, `subphylum`, `superclass`, `class`, `subclass`, `superorder`, `order`, `suborder`, `superfamily`, `family`, `subfamily`, `tribe`, `genus`, `subgenus`, `species group`, `species`, `varietas`, `subspecies`, `serogroup`, `biotype`, `strain`, `serotype`, `genotype`, `forma`, `forma specialis`, `isolate`, `clade`, `no rank`, `subkingdom`, `section`, `FILE`, `ID`, `leaf`|`genus,species group,subgenus`|The ranks or tax ids for which the taxonomy tree is supposed to be refined.|`kmerindexbloom`, `dendrolatex`, `intersectcount`, `intersectcsv`|
