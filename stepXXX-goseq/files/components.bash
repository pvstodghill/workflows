export RESULTS=results
HOWTO="./howto -f ./files/howto.yaml -q"

# ------------------------------------------------------------------------

function reset_results
{
    rm -rf $RESULTS
    mkdir -p $RESULTS
}

# ------------------------------------------------------------------------

function notice_deseq2_regions
{
    echo 1>&2 "# ${FUNCNAME[0]} $*"
    deseq2_dir="$1" ; shift 1
    regions_gff=$deseq2_dir/regions.gff
    regions_tag=locus_tag
    cp $regions_gff $RESULTS/genes.gff
    cat $RESULTS/genes.gff \
	| ./files/gff-to-go-genes \
	      $regions_tag \
	      > ${RESULTS}/genes.txt
    cut -f1 < ${RESULTS}/genes.txt  > ${RESULTS}/genes_assayed.txt
    cut -f2 < ${RESULTS}/genes.txt  > ${RESULTS}/genes_lengths.txt
}
    
# ------------------------------------------------------------------------

function notice_go_download
{
    echo 1>&2 "# ${FUNCNAME[0]} $*"
    go_dir="$1" ; shift 1
    ./files/make-go-funcs.pl \
	$RESULTS/genes.gff \
	$go_dir \
	> ${RESULTS}/go-funcs.txt
}


# ------------------------------------------------------------------------

function notice_deseq2_regulon
{
    echo 1>&2 "# ${FUNCNAME[0]} $*"
    regulon_gff="$1" ; shift 1
    cat $regulon_gff \
	| ./files/make-regulon.pl $RESULTS/genes.gff /dev/stdin \
					> ${RESULTS}/regulon.txt
    cat $regulon_gff \
	| fgrep '; colour 3; ' \
	| ./files/make-regulon.pl $RESULTS/genes.gff /dev/stdin \
					> ${RESULTS}/regulon_up.txt
    cat $regulon_gff \
	| fgrep '; colour 2; ' \
	| ./files/make-regulon.pl $RESULTS/genes.gff /dev/stdin \
					> ${RESULTS}/regulon_down.txt

}

# ------------------------------------------------------------------------

function run_goseq
{
    echo 1>&2 "# ${FUNCNAME[0]} $*"
    label=$1 ; shift 1

    for x in "" _up _down ; do
	$HOWTO Rscript ./files/run-goseq.R "$@" \
	       ${RESULTS}/genes_assayed.txt ${RESULTS}/genes_lengths.txt \
	       ${RESULTS}/regulon${x}.txt ${RESULTS}/go-funcs.txt \
	       > ${RESULTS}/results${x}.${label}.tsv
    done
}


