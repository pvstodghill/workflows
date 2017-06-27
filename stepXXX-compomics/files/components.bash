export HOWTO="../howto -f ./files/howto.yaml -q"

function _ensure {
    for f in "$@" ; do
	if [ ! -e "$f" ] ; then
	    echo 1>&2 "Not found: $f"
	    exit 1
	fi
    done
}

function _preclude {
    for f in "$@" ; do
	if [ -e "$f" ] ; then
	    echo 1>&2 "Found: $f"
	    exit 1
	fi
    done
}

# ------------------------------------------------------------------------
# User modifiable configuration variables
# ------------------------------------------------------------------------

# Parameters for SearchCLI, including search engines
#SEARCHCLI_ARGS="-comet 1 -threads 1"

# Parameters for PeptideShakerCLI
#PEPTIDESHAKERCLI_ARGS="-threads 1"

if [ -z ${RESULTS} ] ; then
    RESULTS=./results
fi

# something like this is a good idea.
#export JVM_FLAGS = -d64 -Xms512m -Xmx16g

# ------------------------------------------------------------------------
# Reset the work directory and output files
# ------------------------------------------------------------------------

function reset_results {
    echo + $FUNCNAME "$@"
    ${HOWTO} -p SearchCLI PeptideShakerCLI getorf msconvert
    rm -rf ${RESULTS}
    mkdir -p ${RESULTS}
}

# ------------------------------------------------------------------------
# Routines for adding protein databases
# ------------------------------------------------------------------------

# make_6frame_protein_database MIN_ORF_SIZE_AA .../genome/*.fna
# - 6-frame translation

function make_6frame_protein_database {
    echo + $FUNCNAME "$@"
    _preclude ${RESULTS}/db.fasta
    MIN_ORF_SIZE_AA=$1
    if [ -z "$MIN_ORF_SIZE_AA" ] ; then
	MIN_ORF_SIZE_AA=10
    fi
    MIN_ORF_SIZE_NC=$[ ${MIN_ORF_SIZE_AA} * 3 ]
    ./files/pvs-cat "$@" >> ${RESULTS}/genome.fna
    make-morfs -s ${MIN_ORF_SIZE_NC} \
	  -f ${RESULTS}/db.fasta -g ${RESULTS}/db.gff \
	  ${RESULTS}/genome.fna
}

# add_raw_protein_database .../genome/*.faa
# - raw protein database

function add_raw_protein_database {
    echo + $FUNCNAME "$@"
    _preclude ${RESULTS}/db.fasta
    # ./files/pvs-cat "$@" \
    # 	| sed -r -e 's/(>[^ ]+) .*/\1/' \
    # 	      >> ${RESULTS}/db.fasta
     ./files/pvs-cat "$@" \
     	| sed -r -e "s/'/\&apos;/g" \
     	      >> ${RESULTS}/db.fasta
}

# make_ncbi_protein_database .../genome/refseq.faa .../genome/refseq*.gff
# make_ncbi_protein_database .../genome/genbank.faa .../genome/genbank*.gff
# - make a protein database with human readable gene names and product descriptions

function make_ncbi_protein_database {
    echo + $FUNCNAME "$@"
    _preclude ${RESULTS}/db.fasta
    ./files/make-componics-fasta "$@" > ${RESULTS}/db.fasta
}

# ------------------------------------------------------------------------
# Add decoys to a protein databse
# ------------------------------------------------------------------------

function add_decoys {
    echo + $FUNCNAME "$@"
    _ensure ${RESULTS}/db.fasta
    ${HOWTO} FastaCLI -in ${RESULTS}/db.fasta -decoy
}

# ------------------------------------------------------------------------
# Create the parameters file for SearchGUI and PeptideShaker
# ------------------------------------------------------------------------

function make_params {
    echo + $FUNCNAME "$@"
    _ensure ${RESULTS}/db_concatenated_target_decoy.fasta
    _preclude ${RESULTS}/params.par
    ${HOWTO} IdentificationParametersCLI  \
	  -out ${RESULTS}/params.par \
	  -db ${RESULTS}/db_concatenated_target_decoy.fasta \
	  "$@"
}

# ------------------------------------------------------------------------
# Prepare the spectra files
# ------------------------------------------------------------------------

NUM_SPECTRA=0
declare -a SPECTRA

function _add_spectra {
    IN_PATH="$1" ; shift 1
    IN_EXT="$1" ;  shift 1
    OUT_DIR="$1" ;  shift 1
    OUT_BASE="$1" ;  shift 1
    COPY_PATH="$OUT_DIR/$OUT_BASE.$IN_EXT"
    (
	PS4='++ '
	set -x
	cp "$IN_PATH" "$COPY_PATH"
	${HOWTO} msconvert --mgf -o "$OUT_DIR" --outfile "$OUT_BASE.mgf" "$COPY_PATH"
    )
    
}

function add_spectra {
    echo + $FUNCNAME "$@"
    if [ -z "$1" -o "$2" ] ; then
	echo 1>&2 "Call add_spectra with exactly spectra file"
	exit 1
    fi
    i=$NUM_SPECTRA
    NUM_SPECTRA=$[$NUM_SPECTRA + 1]
    mkdir -p ${RESULTS}/spectrum
    OUT_DIR=${RESULTS}/spectrum
    OUT_BASE=spectra${i}
    OUT_PATH=$OUT_DIR/$OUT_BASE.mgf
    case "$1" in
	*.mzML) _add_spectra "$1" mzML $OUT_DIR $OUT_BASE ;;
	*.mzML.gz) _add_spectra "$1" mzML.gz $OUT_DIR $OUT_BASE ;;
	*.mgf)
	    (
		PS4='++ '
		set -x
		cp "$1" $OUT_PATH
	    )
	    ;;
	*)
	    echo 1>&2 "Do not know how to handle $1"
	    exit 1
    esac
    SPECTRA[$i]=$OUT_PATH
}

# ------------------------------------------------------------------------
# Run SearchGUI and each of the spectra files
# ------------------------------------------------------------------------

function run_searchcli {
    echo + $FUNCNAME "$@"
    _ensure ${RESULTS}/params.par
    if [ ${#SPECTRA[*]} = 0 ] ; then
	echo 1>&2 "No spectra"
	exit 1
    fi
    for i in `seq 0 $[ ${#SPECTRA[*]} - 1]` ; do
	rm -rf ${RESULTS}/searchruns/${i} ${RESULTS}/searchruns/${i}.*
	mkdir -p ${RESULTS}/searchruns/${i}
	${HOWTO} SearchCLI \
	      -spectrum_files ${SPECTRA[$i]} \
	      -output_folder ${RESULTS}/searchruns/${i} \
	      -id_params ${RESULTS}/params.par \
	      "$@"
	mv ${RESULTS}/searchruns/${i}/searchgui_out.zip \
	   ${RESULTS}/searchruns/${i}.searchgui_out.zip
    done
}


# ------------------------------------------------------------------------
# Run PeptideShaker
# ------------------------------------------------------------------------

function run_peptideshakercli {
    echo + $FUNCNAME "$@"
    _ensure ${RESULTS}/searchruns
    rm -f ${RESULTS}/peptideshaker_out.cpsx
    ${HOWTO} PeptideShakerCLI \
		-experiment experiment \
		-sample sample \
		-replicate 1 \
		-spectrum_files ${RESULTS}/spectrum \
		-identification_files ${RESULTS}/searchruns \
		-id_params ${RESULTS}/params.par \
		-out ${RESULTS}/peptideshaker_out.cpsx \
		"$@"
	rm -rf derby.log resources
}

# ------------------------------------------------------------------------
# Generate the mzIdentML file
# ------------------------------------------------------------------------

function make_mzid_usage {
    echo 1>&2 "Usage: $0 ..."
    echo 1>&2 "-A STR - org_address [contact_address]"
    echo 1>&2 "-E STR - org_email [contact_email]"
    echo 1>&2 "-N STR - org_name [contact_first contact_last]"
    echo 1>&2 "-a STR - contact_address"
    echo 1>&2 "-e STR - contact_email"
    echo 1>&2 "-f STR - contact_first"
    echo 1>&2 "-h - this message"
    echo 1>&2 "-l STR - contact_last"
    exit "$@"
}    


function make_mzid {
    echo + $FUNCNAME "$@"
    _ensure ${RESULTS}/peptideshaker_out.cpsx
    #
    local OPTIND
    local opt_A
    local opt_E
    local opt_N
    local opt_a="missing_address"
    local opt_e="missing_email"
    local opt_f="missing_first_name"
    local opt_h
    local opt_l="missing_last_name"
    while getopts 'A:E:N:a:e:f:hl:' opt ; do
	case "$opt" in
	    A) opt_A="$OPTARG" ;;
	    E) opt_E="$OPTARG" ;;
	    N) opt_N="$OPTARG" ;;
	    a) opt_a="$OPTARG" ;;
	    e) opt_e="$OPTARG" ;;
	    f) opt_f="$OPTARG" ;;
	    h) opt_h=1 ;;
	    l) opt_l="$OPTARG" ;;

	    \?) make_mzid_usage 1 ;;
	    *) echo "Can't happen" ; exit 1 ;;
	esac
    done
    shift $((OPTIND-1))
    if [ "$opt_h" ] ; then
	make_mzid_usage
    fi
    if [ -z "$opt_A" ] ; then
	opt_A="$opt_a"
    fi
    if [ -z "$opt_E" ] ; then
	opt_E="$opt_e"
    fi
    if [ -z "$opt_N" ] ; then
	opt_N="$opt_f $opt_l"
    fi
    #
    rm -f ${RESULTS}/results.mzid
    ${HOWTO} MzidCLI \
	     -in ${RESULTS}/peptideshaker_out.cpsx \
	     -contact_first_name "$opt_f" \
	     -contact_last_name "$opt_l" \
	     -contact_email "$opt_e" \
	     -contact_address "$opt_a" \
	     -organization_name "$opt_N" \
	     -organization_email "$opt_E" \
	     -organization_address "$opt_A" \
	     -output_file ${RESULTS}/results.mzid "$@"
    rm -rf derby.log resources
}

# ------------------------------------------------------------------------
# End of file
# ------------------------------------------------------------------------

