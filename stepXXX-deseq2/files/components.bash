export RESULTS=results
HOWTO="../howto -f ../howto.yaml -q"

# ------------------------------------------------------------------------

function reset_results
{
    rm -rf $RESULTS
    mkdir -p $RESULTS
}

# ------------------------------------------------------------------------

function add_regions
{
    echo 1>&2 "# add_regions $*"
    local gff_file; if [ "$1" ] ; then gff_file=$1 ; shift 1 ; fi
    local feature_type; if [ "$1" ] ; then feature_type=$1 ; shift 1 ; fi
    if [ "$1" ] ; then
	echo 1>&2 "add_regions: multiple features not implemented"
	exit 1
    elif [ -z "$feature_type" ] ; then
	cat "$gff_file" >> $RESULTS/regions.gff
    else
	cat "$gff_file" | fgrep "	$feature_type	" >> $RESULTS/regions.gff

    fi
}

# ------------------------------------------------------------------------

function add_profile
{
    echo 1>&2 "# add_profiles $*"
    local condition; if [ "$1" ] ; then condition=$1 ; shift 1 ; fi
    local sample; if [ "$1" ] ; then sample=$1 ; shift 1 ; fi
    local profile; if [ "$1" ] ; then profile=$1 ; shift 1 ; fi
    echo $sample:$profile >> ${RESULTS}/_counts_args.txt
    echo $condition:$sample >> ${RESULTS}/_prep_args.txt
}

# ------------------------------------------------------------------------

function run_deseq2
{
    echo 1>&2 "# run_deseq2 $*"

    local counts_args="`cat ${RESULTS}/_counts_args.txt`"
    local prep_args="`cat ${RESULTS}/_prep_args.txt`"

    cat $RESULTS/regions.gff \
	| ./files/make-deseq-counts $counts_args \
					  > $RESULTS/counts.txt

    ./files/prep-deseq -q -x -d $RESULTS -c $RESULTS/counts.txt \
			     "$@" $prep_args

    $HOWTO Rscript ./files/run-deseq2 $RESULTS/params.R

    cat $RESULTS/output-extended.txt \
	| ./files/deseq-output2results > $RESULTS/results.txt

}

# ------------------------------------------------------------------------

function make_gff_results
{
    echo 1>&2 "# make_gff_results $*"
    local chrom; if [ "$1" ] ; then chrom=$1 ; shift 1 ; fi
    local cutoff; if [ "$1" ] ; then cutoff=$1 ; shift 1 ; fi
    cat $RESULTS/output.txt \
	| ./files/deseq-output2gff ${chrom} $cutoff \
					 > $RESULTS/${chrom}_results_${cutoff}.gff
    cat $RESULTS/${chrom}_results_${cutoff}.gff \
	| egrep '; colour [23];' > $RESULTS/${chrom}_changed_${cutoff}.gff
}

