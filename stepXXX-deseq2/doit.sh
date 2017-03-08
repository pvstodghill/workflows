#! /bin/bash

set -e

export PATH=/usr/local/bin:/usr/bin:/bin

cd `dirname $0`

. ./files/components.bash

# ------------------------------------------------------------------------

echo 1>&2 "Edit $0"
exit 1

reset_results

ROOT=/some/path
CHROMOSOME=chromosome

add_regions ${CHROMOSOME}.gff gene

add_profile before before-1 $ROOT/${CHROMOSOME}_before-1.unique.sinister.profile
add_profile before before-2 $ROOT/${CHROMOSOME}_before-2.unique.sinister.profile
add_profile before before-3 $ROOT/${CHROMOSOME}_before-3.unique.sinister.profile
add_profile after after-1 $ROOT/${CHROMOSOME}_after-1.unique.sinister.profile
add_profile after after-2 $ROOT/${CHROMOSOME}_after-2.unique.sinister.profile
add_profile after after-3 $ROOT/${CHROMOSOME}_after-3.unique.sinister.profile

run_deseq2

make_gff_results ${CHROMOSOME} 2.0
make_gff_results ${CHROMOSOME} 1.5
