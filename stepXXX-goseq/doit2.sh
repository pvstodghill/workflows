#! /bin/bash

set -e

export PATH=/usr/local/bin:/usr/bin:/bin

cd `dirname $0`
. ./files/components.bash

# ------------------------------------------------------------------------

echo 1>&2 "Edit $0"
exit 1

# run_goseq_once FIXME/regulon.txt strict  -q 0.01
# run_goseq_once FIXME/regulon.txt relaxed -q 0.05
# run_goseq_once FIXME/regulon.txt all     -q 1

run_goseq strict  -q 0.01
run_goseq relaxed -q 0.05
run_goseq all     -q 1
