#!/bin/sh
set -e

scriptdir=$(dirname "$0")
cd $scriptdir/../../../../..

basedir=$(pwd)

commondir=${basedir}/data/common
kudir=${basedir}/krakenuniq

export JELLYFISH_BIN=${kudir}/jellyfish-install/bin/jellyfish

mkdir -p ./viral_db
mkdir -p ./viral_db/library
mkdir -p ./viral_db/taxonomy

cp ${commondir}/nodes.dmp ./viral_db/taxonomy
cp ${commondir}/names.dmp ./viral_db/taxonomy

# KU only accepts unzipped files ending with fa, fna and so on
cp ./data/projects/viral/fasta/*.fa ./viral_db/library
cp ./data/projects/viral/csv/viral_ku.map ./viral_db/library

${kudir}/krakenuniq-build --kmer-len 31 --db ./viral_db
# Only one thread here so that results are in deterministic order for later comparison
${kudir}/krakenuniq --threads 1 --db ./viral_db --output ./core/target/test.fasta.out ./core/src/test/resources/projects/viral/test.fasta.gz

mkdir -p ./viral16_db
mkdir -p ./viral16_db/library
mkdir -p ./viral16_db/taxonomy

cp ${commondir}/nodes.dmp ./viral16_db/taxonomy
cp ${commondir}/names.dmp ./viral16_db/taxonomy

# KU only accepts unzipped files ending with fa, fna and so on
cp ./data/projects/viral/fasta/*.fa ./viral16_db/library
cp ./data/projects/viral/csv/viral_ku.map ./viral16_db/library

${kudir}/krakenuniq-build --kmer-len 16 --db ./viral16_db
# Only one thread here so that results are in deterministic order for later comparison
${kudir}/krakenuniq --threads 1 --db ./viral16_db --output ./core/target/test16.fasta.out ./core/src/test/resources/projects/viral/test.fasta.gz
