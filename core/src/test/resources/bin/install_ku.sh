#!/bin/sh
set -e

scriptdir=$(dirname "$0")
cd $scriptdir/../../../../..

basedir=$(pwd)
kudir=${basedir}/krakenuniq

# Check whether a usable local KrakenUniq installation already exists.
# We deliberately do more than just test for the krakenuniq directory: a
# leftover, partial or broken build would still leave the directory in place.
# Instead we verify that every executable actually used by
# make_n_apply_viral_dbs.sh is present, executable and genuinely runnable
# (i.e. it returns its version without error).
ku_installed() {
    krakenuniq="${kudir}/krakenuniq"
    krakenuniq_build="${kudir}/krakenuniq-build"
    jellyfish="${kudir}/jellyfish-install/bin/jellyfish"

    for bin in "$krakenuniq" "$krakenuniq_build" "$jellyfish"; do
        [ -x "$bin" ] || return 1
    done

    # Actually invoke the binaries so we detect a broken/incomplete build.
    "$krakenuniq" --version >/dev/null 2>&1 || return 1
    "$krakenuniq_build" --version >/dev/null 2>&1 || return 1
    "$jellyfish" --version >/dev/null 2>&1 || return 1

    return 0
}

# Set FORCE_KU_INSTALL=1 to always (re-)install, even if a usable installation exists.
if [ "${FORCE_KU_INSTALL:-0}" = "1" ]; then
    echo "FORCE_KU_INSTALL=1 set - (re-)installing KrakenUniq into ${kudir} ..."
elif ku_installed; then
    echo "KrakenUniq already installed and runnable in ${kudir} - skipping installation."
    exit 0
else
    echo "No usable KrakenUniq installation found - installing into ${kudir} ..."
fi

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
