#!/bin/sh
scriptdir=$(dirname "$0")

java -Xmx32g -jar $scriptdir/../lib/genestrip.jar -d $scriptdir/../data "$@" 