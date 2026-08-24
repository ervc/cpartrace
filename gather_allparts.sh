#!/usr/bin/env bash
__usage="
Read allparts_rank* files and create one compiled allparts.txt
Usage: $0 outputdir
"

if [[ "$#" -ne 1 ]]; then
    echo "$__usage"
    exit 1
fi
if [[ "#1" == "-h" ]]; then
    echo "$__usage"
    exit 0
fi

# list allparts from directory
OUTDIR="$1"
LASTRANK=$(find "$OUTDIR"/allparts_rank*.txt | xargs -- basename -a | sed s/[^0-9]//g | sort -n | tail -1)
echo $LASTRANK
