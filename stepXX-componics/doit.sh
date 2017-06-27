#! /bin/bash

set -e

#export PATH=/usr/local/bin:/usr/bin:/bin

cd `dirname $0`

export JVM_FLAGS="-d64 -Xms512m -Xmx16g"

. files/components.bash

reset_results

# v--- pick one ---v

# MIN_ORF_SIZE_AA=10
# make_protein_database $MIN_ORF_SIZE_AA .../genome/*.fna

# add_raw_protein_database .../proteome.faa

make_ncbi_protein_database .../genome/proteome.faa .../genome/annotation*.gff

# ^--- pick one ---^

add_decoys

FIXED_MODS="Methylthio of C"
VARIABLE_MODS="Oxidation of M"
VARIABLE_MODS+=", Deamidation of N"
VARIABLE_MODS+=", Deamidation of Q"
VARIABLE_MODS+=", Formylation of protein N-term"

# -comet_enzyme_type 1 == semi-specific

rm -f ./results/params.par
make_params \
    -mc 1 \
    -frag_tol 0.6 \
    -fixed_mods "${FIXED_MODS}" \
    -variable_mods "${VARIABLE_MODS}" \
    -comet_min_peaks 1 \
    -comet_enzyme_type 1 \
    -comet_min_prec_mass 350.0 \
    -comet_max_prec_mass 5000.0 \
    #

add_spectra .../spectra/*.mgf

SEARCHCLI_ARGS=
#SEARCHCLI_ARGS+=" -threads 1"
SEARCHCLI_ARGS+=" -comet 1"

run_searchcli $SEARCHCLI_ARGS

PEPTIDESHAKERCLI_ARGS=
#PEPTIDESHAKERCLI_ARGS+=" -threads 1"

run_peptideshakercli $PEPTIDESHAKERCLI_ARGS

make_mzid -f Paul -l Stodghill \
	  -e paul.stodghill@ars.usda.gov \
	  -a "Robert Holley Center, 538 Tower Road, Ithaca, NY 14853"
