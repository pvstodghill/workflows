HOWTO="../howto -f ../howto.yaml -q"

# ------------------------------------------------------------------------
# Configuration variables
# ------------------------------------------------------------------------

### proteins.faa used to produce .pep.xml files
if [ -z "${DB_FASTA}" ] ; then
    echo 1>&2 "DB_FASTA not set."
    exit 1
fi

### cannonical annotation.gff for genome
if [ -z "${ANNOT_GFF}" ] ; then
    echo 1>&2 "ANNOT_GFF not set."
    exit 1
fi

### directory in which to write intermediate files and results
: ${RESULTS:=results}
export RESULTS

### PSM cutoff
: ${PEP_XML_QVALUE:=0.01}

### adjust error rate for < this number of false positives
: ${FALSE_POSITIVES_THRESHOLD:=1}

### min number of peptides needed to support an ORF
: ${ORF_SUPPORT:=2}

### min number of sample ORF's needed to support an ORF
: ${SAMPLE_SUPPORT:=2}


# ------------------------------------------------------------------------
# fixme
# ------------------------------------------------------------------------

function reset_results {
    echo + ${FUNCNAME[0]} "$@"

    rm -rf ${RESULTS}
    rm -f results.zip
    mkdir -p ${RESULTS}
}


# ------------------------------------------------------------------------
# fixme
# ------------------------------------------------------------------------

SAMPLES=

function add_sample {
    echo + ${FUNCNAME[0]} "$@"

    I=$1 ; shift 1


    ${HOWTO} perl ./files/pep-xml-to-psm.pl  -q ${PEP_XML_QVALUE} -l ${I} "$@"  > ${RESULTS}/${I}_raw.psms.txt

    cat ${RESULTS}/${I}_raw.psms.txt \
	| ${HOWTO} perl ./files/psm-errors.pl  \
		  -c ${FALSE_POSITIVES_THRESHOLD} -i q_value \
		  > ${RESULTS}/${I}_tight.psms.txt


    cat ${RESULTS}/${I}_tight.psms.txt \
	| ${HOWTO} perl ./files/psm-map.pl  -p 1 ${DB_FASTA} \
		  > ${RESULTS}/${I}_mapped.psms.txt


    cat ${RESULTS}/${I}_mapped.psms.txt \
	| ${HOWTO} perl ./files/psm-to-peptide-gff.pl -p  ${DB_FASTA} \
		  > ${RESULTS}/peptides_${I}.gff


    cat ${RESULTS}/${I}_mapped.psms.txt \
	| ${HOWTO} perl ./files/psm-to-protein-gff.pl  \
		  -o -p ${ORF_SUPPORT} ${DB_FASTA}  ${ANNOT_GFF} \
		  > ${RESULTS}/orfs_${I}.gff


    cat ${RESULTS}/peptides_${I}.gff | cat \
	| ${HOWTO} perl ./files/split-gff.pl -n -d ${RESULTS} peptides_${I}


    cat ${RESULTS}/orfs_${I}.gff | cat \
	| ${HOWTO} perl ./files/split-gff.pl -n -d ${RESULTS} orfs_${I}

    SAMPLES+=" ${I}"
}

# ------------------------------------------------------------------------
# fixme
# ------------------------------------------------------------------------

function make_venn_diagram {
    echo + ${FUNCNAME[0]} "$@"

    ORFS_FILES=
    for I in $SAMPLES ; do
	ORFS_FILES+=" ${RESULTS}/orfs_${I}.gff"
    done

    ${HOWTO} perl ./files/make-venn.pl -o ${RESULTS}/venn.png $ORFS_FILES
}

# ------------------------------------------------------------------------
# hv_orfs.txt - list of ORF's supported by >=${SAMPLE_SUPPORT} samples
# ------------------------------------------------------------------------

function make_hv_orfs {
    echo + ${FUNCNAME[0]} "$@"

    ORFS_FILES=
    for I in $SAMPLES ; do
	ORFS_FILES+=" ${RESULTS}/orfs_${I}.gff"
    done

    cat $ORFS_FILES \
	| egrep -v '^#' \
	| sed -e 's/.*orf_id=//' -e 's/;.*//'  \
	| sort | uniq -c \
	| gawk 'OFS=""; {if ($1>='${SAMPLE_SUPPORT}') { print "protein:",$2; }}' \
	       > ${RESULTS}/hv_orfs.txt

}

# ------------------------------------------------------------------------
# hv_mapped.psms.txt - set of high-value PSM's
# ------------------------------------------------------------------------

function make_hv_mapped_psms {
    
    echo + ${FUNCNAME[0]} "$@"

    ONE_MAPPED_FILE=
    MAPPED_FILES=
    for I in $SAMPLES ; do
	ONE_MAPPED_FILE=${RESULTS}/${I}_mapped.psms.txt
	MAPPED_FILES+=" ${RESULTS}/${I}_mapped.psms.txt"
    done


    (
	head -n1 ${ONE_MAPPED_FILE}
	cat $MAPPED_FILES \
	    | fgrep -f ${RESULTS}/hv_orfs.txt
    )  > ${RESULTS}/hv_mapped.psms.txt
}

# ------------------------------------------------------------------------
# *peptides_hv.gff - set of high-value (mapped) peptides
# ------------------------------------------------------------------------

function make_peptide_hv_gff {
    echo + ${FUNCNAME[0]} "$@"

    cat ${RESULTS}/hv_mapped.psms.txt \
	| ${HOWTO} perl ./files/psm-to-peptide-gff.pl -p  ${DB_FASTA} \
		  > ${RESULTS}/peptides_hv.gff

    cat ${RESULTS}/peptides_hv.gff \
	| ${HOWTO} perl ./files/split-gff.pl -n -d results peptides_hv
}

# ------------------------------------------------------------------------
# *regions_hv.gff - set of regions supported by high-value (mapped) peptides
# ------------------------------------------------------------------------

function make_regions_hv_gff {
    echo + ${FUNCNAME[0]} "$@"

    cat ${RESULTS}/hv_mapped.psms.txt \
	| ${HOWTO} perl ./files/psm-to-protein-gff.pl  \
		  -r -p 2 ${DB_FASTA}  ${ANNOT_GFF} \
		  > ${RESULTS}/regions_hv.gff

    cat ${RESULTS}/regions_hv.gff | cat \
	| ${HOWTO} perl ./files/split-gff.pl -n -d results regions_hv
}

# ------------------------------------------------------------------------
# *overlaps_hv.gff - ???
# *opposites_hv.gff - ???
# ------------------------------------------------------------------------

function make_other_hv_gff {
    echo + ${FUNCNAME[0]} "$@"

    cat ${RESULTS}/*_regions_hv.gff \
	| ${HOWTO} perl ./files/flag-overlapping.pl  -a overlaps_regions -o \
	| ${HOWTO} perl ./files/split-gff.pl -n -d results overlaps_hv

    cat ${RESULTS}/*_overlaps_hv.gff \
	| fgrep overlaps_regions_opposite \
	| ${HOWTO} perl ./files/split-gff.pl -n -d results opposites_hv

}

# ------------------------------------------------------------------------
# fixme
# ------------------------------------------------------------------------

function process_samples {
    make_hv_orfs
    make_hv_mapped_psms
    make_peptide_hv_gff
    make_regions_hv_gff
    make_other_hv_gff
}



# ------------------------------------------------------------------------
# Print some stats
# ------------------------------------------------------------------------

function print_some_stats {
    I="$1"
    TAG="$2"
    OBJS="$3"
    echo "= ${TAG} ="
    echo -n "peptides: "
    cat ${RESULTS}/*_peptides_${I}.gff | egrep -v '^#' | wc -l
    echo -n "- max p-value: "
    cat ${RESULTS}/*_peptides_${I}.gff | sed -e 's/.*PEP=//' -e 's/;.*//' -e 's/"//g' | sort -rg | head -n1
    echo -n "- max q-value: "
    cat ${RESULTS}/*_peptides_${I}.gff | sed -e 's/.*qvalue=//' -e 's/;.*//' -e 's/"//g' | sort -rg | head -n1
    echo -n "supported ${OBJS}: "
    cat ${RESULTS}/*_${OBJS}_${I}.gff | egrep -v '^#' | wc -l
    echo -n "- max p-value: "
    cat ${RESULTS}/*_${OBJS}_${I}.gff | sed -e 's/.*pvalue=//' -e 's/;.*//' -e 's/"//g' | sort -rg | head -n1
    echo ""
}

    
function print_stats {
    echo + ${FUNCNAME[0]} "$@"

    echo ""
    for I in $SAMPLES ; do
	print_some_stats "${I}" "${I}" orfs
    done
    print_some_stats hv High-value regions

    
    echo -n "... without CDSs: "
    cat ${RESULTS}/*_regions_hv.gff | fgrep missing_cds | wc -l
    echo -n "... with CDSs requiring extension:"
    cat ${RESULTS}/*_regions_hv.gff | fgrep reannotate_cds | wc -l
    echo -n "... with pseudogenes: "
    cat ${RESULTS}/*_regions_hv.gff | fgrep pseudo | wc -l
    echo ""

    echo -n "... that overlap, either strand: "
    cat ${RESULTS}/*_overlaps_hv.gff | egrep -v '^#' | wc -l
    echo -n "... that overlap, opposite strands: "
    cat ${RESULTS}/*_opposites_hv.gff | egrep -v '^#' | wc -l
    echo ""

}
