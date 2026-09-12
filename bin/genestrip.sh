#!/bin/sh
set -e

scriptdir=$(dirname "$0")

# Heap. The largest database in ft-db-exp2 is nocardia, whose fill needs about 37.4 GiB live:
# 18.0 GiB of long[] buckets (RadixKMerStore:294) + ~4.5 GiB of values + 14.9 GiB for the fill
# Bloom filter at fillBloomFilterFpp=1e-11, all three live at once because FillDBGoal:118
# passes the filter as a constructor argument. 32g was below that and FillDBGoal died with a
# clean "OutOfMemoryError: Java heap space". 48g clears it and leaves a 64 GB machine ~16 GB
# for the OS and the page cache, which Genestrip leans on heavily while reading FASTA.
# Override per run with GS_XMX if a machine has more or less.
: "${GS_XMX:=56g}"

java -Xmx"$GS_XMX" -jar $scriptdir/../lib/genestrip.jar -d $scriptdir/../data "$@" 