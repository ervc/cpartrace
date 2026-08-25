#!/usr/bin/env bash
__usage="Read allparts_rank* files and create one compiled allparts.txt
Usage: $0 outputdir [dir2 dir3 ... dirN]
"

if [[ "$#" -eq 0 ]]; then
    echo "$__usage"
    exit 1
fi
if [[ "#1" == "-h" ]]; then
    echo "$__usage"
    exit 0
fi

gather_files () {
    # list allparts from directory
    OUTDIR="$1"
    echo "----------------------------------------"
    echo "Working on $OUTDIR"
    OUTFILE="$OUTDIR"/allparts.txt
    echo "t x y z vx vy vz status" > "$OUTFILE"
    # Note that order is NOT preserved here but that does not matter
    FILES=$(find "$OUTDIR"/allparts_rank*.txt)
    for f in ${FILES[@]}; do
        # trim first line from file
        tail -n +2 "$f" >> "$OUTFILE"
        rm $f
    done
    echo "Saved to $OUTFILE"
    echo "----------------------------------------"
}

for DIR in "$@"; do
    gather_files "$DIR"
done
