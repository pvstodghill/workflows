#! /bin/bash -x

cd `dirname $0`

echo 1>&2 "Edit $0"
exit 1

set -e

rm -f README.geneontology.org
rm -f README.uniprot-idmapping
rm -f go_*-termdb-tables.tar.gz
rm -f idmapping_selected.tab.gz
rm -f term.txt
rm -f term2term.txt


# ------------------------------------------------------------------------

# Script to download most recent GO annotation and DC3000 genome

# Get the current GO annotation from UniProt
IDMAPPING_URL=ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping
wget ${IDMAPPING_URL}/idmapping_selected.tab.gz
wget ${IDMAPPING_URL}/README
mv README README.uniprot-idmapping

GO_URL=http://archive.geneontology.org/latest-full
VERSION=monthly
wget ${GO_URL}/go_${VERSION}-termdb-tables.tar.gz
wget ${GO_URL}/README
mv README README.geneontology.org

# Extract the important tables from the GO tarballs
tar zxf go_${VERSION}-termdb-tables.tar.gz \
    go_${VERSION}-termdb-tables/term.txt \
    go_${VERSION}-termdb-tables/term2term.txt
mv go_${VERSION}-termdb-tables/* .
rmdir go_${VERSION}-termdb-tables


