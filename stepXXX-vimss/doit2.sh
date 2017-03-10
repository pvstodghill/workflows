#! /bin/bash -x

# ------------------------------------------------------------------------
# Convert to GFF
# ------------------------------------------------------------------------

./files/make-vimss.pl results/genes.txt results/named.txt > results/vimss.gff

./files/split-gff -d results vimss < results/vimss.gff
