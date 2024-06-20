#!/bin/sh
set -e

scriptdir=$(dirname "$0")

java -Xmx32g -jar $scriptdir/../lib/genestrip.jar -d $scriptdir/../data "$@" 