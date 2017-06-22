#! /bin/bash

set -e

export PATH=/usr/local/bin:/usr/bin:/bin

cd `dirname $0`

DB_FASTA=extended.fasta
ANNOT_GFF=genbank.gff

. files/components.bash

reset_results

DIR=/path/to/somewhere

add_sample 1 ${DIR}/pellet.1.pep.xml ${DIR}/super.1.pep.xml
add_sample 2 ${DIR}/pellet.2.pep.xml ${DIR}/super.2.pep.xml
add_sample 3 ${DIR}/pellet.3.pep.xml ${DIR}/super.3.pep.xml

process_samples

make_venn_diagram 

print_stats

set -x
zip results.zip \
    ${RESULTS}/venn.png \
    ${RESULTS}/*_peptides_1.gff \
    ${RESULTS}/*_peptides_2.gff \
    ${RESULTS}/*_peptides_3.gff \
    ${RESULTS}/*_regions_hv.gff
