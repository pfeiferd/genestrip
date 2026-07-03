#!/bin/sh
set -e

scriptdir=$(dirname "$0")
cd $scriptdir/../../../../..

basedir=$(pwd)

# We simply use the latest GitHub version... (currently v1.0.4)
rm -rf krakenuniq
git clone https://github.com/fbreitwieser/krakenuniq.git
# if error uint32_t appears, create a fixed version of uid_mapping.hpp with #include <cstdint> and report-cols.hpp in $basedir/bin/fixes/ku/ folder.
# cp $basedir/bin/fixes/ku/*.hpp krakenuniq/src
cd krakenuniq

# This file must be deleted, it leads to a compile error as it is not C++
rm src/gzstream/version

# For fixes under macOS the GNU version of sed is required (BSD sed uses a
# different -i syntax). Install it via: brew install gnu-sed
# On macOS GNU sed is called "gsed", on Linux plain "sed" already is GNU sed.
if command -v gsed >/dev/null 2>&1; then
    SED=gsed
else
    SED=sed
fi

# Fix missing includes:
$SED -i '1i\#include <stdint.h>' src/uid_mapping.hpp
$SED -i '12i\#include <stdint.h>' src/report-cols.hpp

# Now we can build...
./install_krakenuniq.sh -j .
