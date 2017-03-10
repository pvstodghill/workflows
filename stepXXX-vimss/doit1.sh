#! /bin/bash -x

# ------------------------------------------------------------------------
# VIMSS
# - http://www.microbesonline.org/operons/
# - http://www.microbesonline.org/operons/gnc${SOURCE_ID}.html
# ------------------------------------------------------------------------

SOURCE_ID=223283 # DC3000

# ------------------------------------------------------------------------
# Delete old results
# ------------------------------------------------------------------------

rm -rf results
mkdir -p results

# ------------------------------------------------------------------------
# Download VIMSS predictions
# ------------------------------------------------------------------------

wget -O results/genes.txt \
     "http://www.microbesonline.org/cgi-bin/genomeInfo.cgi?tId=${SOURCE_ID};export=tab"
wget -O results/named.txt \
     "http://www.microbesonline.org/operons/gnc${SOURCE_ID}.named"
wget -O results/scores.txt \
     "http://www.microbesonline.org/operons/gnc${SOURCE_ID}.scores"
