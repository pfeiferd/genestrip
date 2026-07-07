#!/bin/sh
set -e

scriptdir=$(dirname "$0")
cd $scriptdir/../../../../..

basedir=$(pwd)

commondir=${basedir}/data/common
kudir=${basedir}/krakenuniq

export JELLYFISH_BIN=${kudir}/jellyfish-install/bin/jellyfish

# Final deliverables of this script (the two output files are installed into the
# local Maven repo by the gentestdata profile).
outfile=${basedir}/core/target/test.fasta.out
outfile16=${basedir}/core/target/test16.fasta.out

# Check whether the viral databases have already been built and applied.
# We deliberately do more than just test for the *_db directories: a leftover,
# partial or aborted run would still leave those directories in place. Instead we
# require the actually built KrakenUniq database file of each db and both
# non-empty classification outputs to be present.
dbs_built() {
    for f in \
        "${basedir}/viral_db/database.kdb" \
        "${basedir}/viral16_db/database.kdb" \
        "$outfile" \
        "$outfile16"; do
        [ -s "$f" ] || return 1
    done
    return 0
}

# Set FORCE_KU_INSTALL=1 to always (re-)build, even if usable databases exist.
if [ "${FORCE_KU_INSTALL:-0}" = "1" ]; then
    echo "FORCE_KU_INSTALL=1 set - (re-)building viral databases ..."
elif dbs_built; then
    echo "Viral databases and outputs already present - skipping build."
    exit 0
else
    echo "No usable viral databases found - building ..."
fi

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
