#!/bin/sh
scriptdir=$(dirname "$0")

java -Xmx16g -jar $scriptdir/../lib/genestrip.jar -d $scriptdir/../data "$@" 