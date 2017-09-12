ROOT=$(cd .. ; /bin/pwd )
export HOWTO="$ROOT/howto -f $ROOT/howto.yaml -q"

RESULTS=results

# ------------------------------------------------------------------------

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

function add_genome_usage {
    echo 1>&2 "Usage: add_genome NAME file1.gbk file2.gbk ..."
    echo 1>&2 "Usage: add_genome NAME file1.fna file2.fna ..."
    exit $*
}


function add_genome {
    NAME="$1" ; shift 1
    if [ -z "$1" ] ; then
	add_genome_usage 1
    fi
    
    mkdir -p ${RESULTS}/genomes

    case "$1" in
	*.gbk) EXT=gbk ;;
	*.fna) EXT=fasta ;;
	*.gbk.tgz) EXT=gbk ;;
	*.fna.tgz) EXT=fasta ;;
	*) add_genome_usage 1
    esac

    rm -f  "${RESULTS}/genomes/${NAME}.${EXT}"
    (
	for f in "$@" ; do
	    case "$f" in
		*.contig.*.tgz) tar zxOf "$f" ;;
		*.*) cat "$f"
	    esac
	done
    )  >> "${RESULTS}/genomes/${NAME}.${EXT}"
}


function run_progressiveMauve {
    $HOWTO progressiveMauve --output=${RESULTS}/alignment.xmfa ${RESULTS}/genomes/*.*
}

function run_Mauve {
    _ensure ${RESULTS}/alignment.xmfa
    $HOWTO Mauve ${RESULTS}/alignment.xmfa
}
