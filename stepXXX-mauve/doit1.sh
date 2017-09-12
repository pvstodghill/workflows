#! /bin/bash

set -e

export PATH=/usr/local/bin:/usr/bin:/bin

cd `dirname $0`
. ./files/components.bash

# ------------------------------------------------------------------------

set -x

rm -rf results
mkdir results

echo 1>&2 "Edit $0"
exit 1


add_genome psyringae_DC3000 /.../*DC3000*/NC_004578.gbk
add_genome psyringae_1448A /.../*1448A*/NC_005773.gbk
add_genome psyringae_B728a /.../*B728a*/NC_007005.gbk

run_progressiveMauve
