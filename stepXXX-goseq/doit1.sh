#! /bin/bash

set -e

export PATH=/usr/local/bin:/usr/bin:/bin

cd `dirname $0`
. ./files/components.bash

# ------------------------------------------------------------------------

echo 1>&2 "Edit $0"
exit 1

reset_results

notice_deseq2_regions ../stepXXX-deseq2/results/regions.gff

notice_go_download FIXME

notice_deseq2_regulon ../stepXXXX-deseq2.results/NC_XXXXXX_changed_N.gff
