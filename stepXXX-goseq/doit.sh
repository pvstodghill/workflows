#! /bin/bash

set -e

export PATH=/usr/local/bin:/usr/bin:/bin

cd `dirname $0`
. ./files/components.bash

# ------------------------------------------------------------------------

echo 1>&2 "Edit $0"
exit 1

reset_results

notice_deseq2_regions step1-deseq2.results

notice_go_download $MIDDEN/2017/01-31-go

notice_deseq2_regulon step1-deseq2.results/NC_002947_changed_1.5.gff

run_goseq strict  -q 0.01
run_goseq relaxed -q 0.05
