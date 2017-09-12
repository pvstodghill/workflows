export HOWTO="../howto -f ../howto.yaml -q"

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

