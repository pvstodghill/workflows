#! /bin/bash

set -e

export PATH=/usr/local/bin:/usr/bin:/bin

cd `dirname $0`
. ./files/components.bash

# ------------------------------------------------------------------------

echo 1>&2 "Edit $0"
exit 1

