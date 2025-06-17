|Name|Description|
|-|-|
|`pos`|Sort position of entry.|
|`name`|The name associated with the tax id.|
|`rank`|The rank of the tax id.|
|`taxid`|The tax id.|
|`reads`|The number of reads classified with respect to the tax id.|
|`kmers from reads`|The number of *k*-mers from classified reads which are consistent with the read's tax id.|
|`kmers`|*All* matched *k*-mers which are specific to the tax id's genome (according to the database). The *k*-mers do not have to be in an accordingly classified read for this.|
|`unique kmers`|*All* unique *k*-mers, which are specific to the tax id's genome (according to the database). Here, multiple occurrences of the same *k*-mer are only counted once. Genestrip always performs exact counting according to [KrakenUniq's exact counting](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-018-1568-0#Sec8). (Genestrip implements an efficient in-memory storage method for related counts based on a bit vector.)|
|`contigs`|The number of contiguous sequences of *k*-mers that are specific to the tax id's genome.|
|`average contig length`|The average length of contiguous sequences of *k*-mers (as the number of *k*-mers) that are specific to the tax id's genome.|
|`max contig length`|The maximum length of all contiguous sequences of *k*-mers (as the number of *k*-mers) that are specific to the tax id's genome.|
|`reads >=1 kmer`|Reads with at least on *k*-mer of the respective tax id.|
|`reads bps`|The total number of base pairs of reads classified with respect to the tax id.|
|`avg. read length`|The average length of classified reads in base pairs.|
|`db coverage`|The ratio `unique kmers` / u<sub>t</sub>, , where *u<sub>t</sub>* = `db kmers`|
|`exp. unique kmers`|The number of expected unique *k*-mers, which is *u<sub>t</sub> * (1 - (1 - 1/u<sub>t</sub>)*<sup>`kmers`</sup>), where *u<sub>t</sub>* is the number of specific *k*-mers for the tax id in the database.|
|`unique kmers / exp.`|The ratio `unique kmers` / `exp. unique kmers` for the tax id. This should be close to 1 for a consistent match of *k*-mers. ([This paper](https://arxiv.org/pdf/1602.05822.pdf) discusses the corresponding background distribution (of `unique kmers`).)|
|`db kmers`|The number *u<sub>t</sub>* of specific *k*-mers for the tax id in the database.|
|`parent taxid`|The parent tax id.|
|`mean error`|The mean ratio of a classified read's *k*-mers that are not in the database per read's total *k*-mers.|
|`kmer error std. dev.`|The standard deviation of the `mean error`.|
|`mean class error`|The mean ratio of a read's *k*-mers that are not consistent with the read's class per read's total *k*-mers.|
|`class error std. dev.`|The standard deviation of the `mean class error`.|
|`contig len std. dev.`|The standard deviation of the `contig length`.|
|`norm. reads`,`norm. kmers`,`norm. reads bps`,`norm. read >=1 kmer`,`norm. reads kmers`|Normalized value of a respective value type.|
|`acc. reads`,`acc. norm. reads`,`acc. kmers`,`acc. norm. kmers`,`acc. reads bps`,`acc. norm. reads bps`,`acc. read >=1 kmer`,`acc. norm. read >=1 kmer`,`acc. reads kmers`,`acc. norm. reads kmers`|Accumulated value or accumulated normalized valued of a respective value type.|
|`max contig desc.`|The descriptor of a read that holds a contiguous sequence of maximum length (according to the column `max contig length`).|
|`acc. mean error`|The accumulated `mean error`.|
|`acc. error std. dev.`|The standard deviation of the `acc. mean error`.|
|`acc. class mean error`|The accumulated `class mean error`.|
|`acc. class error std. dev.`|The standard deviation of the `acc. class mean error`.|
|`max kmer counts`|The frequencies of the most frequent unique *k*-mers which are specific to the tax id's genome in descending order separated by `;`. This column is experimental and only added when the configuration property `matchWithKMerCounts` is set to `true`. The number of frequencies is determined via `maxKMerResCounts` (see also Section [Configuration parameters](#configuration-parameters)).|
