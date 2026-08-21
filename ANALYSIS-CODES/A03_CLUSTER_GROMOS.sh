#!/usr/bin/env bash

set -euo pipefail
export LC_ALL=C

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd -- "${SCRIPT_DIR}/.." && pwd)"

OUTPUT_DIR="${ROOT_DIR}/03_GROMACS/06_OUTPUT"

ANALYSIS_ROOT="${ROOT_DIR}/04_ANALYSIS"
ANALYSIS_DIR="${ANALYSIS_ROOT}/03_CLUSTER_GROMOS"

mkdir -p "$ANALYSIS_ROOT"
mkdir -p "$ANALYSIS_DIR"

while true; do

    read -rp "RMSD cutoff [nm]: " CUTOFF

    if [[ "$CUTOFF" =~ ^([0-9]+([.][0-9]*)?|[.][0-9]+)$ ]] &&
       awk -v x="$CUTOFF" 'BEGIN { exit !(x > 0) }'
    then
        break
    fi

    echo "ERROR: Please enter a positive number."
    echo "Example: 0.20"

done

while true; do

    read -rp "Analyze every N-th frame: " SKIP

    if [[ "$SKIP" =~ ^[1-9][0-9]*$ ]]; then
        break
    fi

    echo "ERROR: Please enter a positive integer."
    echo "Example: 10"

done

XTC=$(
    find "$OUTPUT_DIR" \
        -maxdepth 1 \
        -type f \
        -name "*.xtc" \
        -printf "%T@ %p\n" \
    | sort -nr \
    | awk '
        NR == 1 {
            $1 = ""
            sub(/^ /, "")
            print
        }
    '
)

PDB=$(
    find "$OUTPUT_DIR" \
        -maxdepth 1 \
        -type f \
        -name "*.pdb" \
        -printf "%T@ %p\n" \
    | sort -nr \
    | awk '
        NR == 1 {
            $1 = ""
            sub(/^ /, "")
            print
        }
    '
)

if [ -z "${XTC:-}" ]; then
    echo "ERROR: XTC file not found in:"
    echo "  $OUTPUT_DIR"
    exit 1
fi

if [ -z "${PDB:-}" ]; then
    echo "ERROR: PDB file not found in:"
    echo "  $OUTPUT_DIR"
    exit 1
fi

echo
echo "========================================"
echo " Clustering Analysis"
echo "========================================"
echo
echo "Method      : GROMOS"
echo "RMSD cutoff : ${CUTOFF} nm"
echo "Frame skip  : ${SKIP}"
echo
echo "========================================"
echo

printf "0\n0\n" | \
gmx cluster \
    -s "$PDB" \
    -f "$XTC" \
    -method gromos \
    -cutoff "$CUTOFF" \
    -skip "$SKIP" \
    -tu ns \
    -g "${ANALYSIS_DIR}/clustering-log.log" \
    -clid "${ANALYSIS_DIR}/cluster-id.xvg" \
    -clndx "${ANALYSIS_DIR}/cluster-index.ndx" \
    -sz "${ANALYSIS_DIR}/cluster-size.xvg" \
    -ntr "${ANALYSIS_DIR}/cluster-trans.xvg" \
    -tr "${ANALYSIS_DIR}/cluster-trans.xpm" \
    -cl "${ANALYSIS_DIR}/representative-structures.pdb" \
    -o "${ANALYSIS_DIR}/rmsd-clust.xpm" \
    -dist "${ANALYSIS_DIR}/rmsd-dist.xvg" \
    -om "${ANALYSIS_DIR}/rmsd-raw.xpm"

awk '
BEGIN {
    print "Time [ns],ID"
}

/^[[:space:]]*#/ {
    next
}

/^[[:space:]]*@/ {
    next
}

NF >= 2 {
    print $1 "," $2
}
' "${ANALYSIS_DIR}/cluster-id.xvg" \
> "${ANALYSIS_DIR}/cluster-id.csv"

rm -f "${ANALYSIS_DIR}"/cluster-*.pdb

awk -v outdir="$ANALYSIS_DIR" '

/^MODEL[[:space:]]+/ {

    model = $2
    outfile = outdir "/cluster-" model ".pdb"

    print > outfile
    next
}

model != "" {

    if (/^ENDMDL/) {

        print > outfile
        print "END" > outfile

        close(outfile)

        model = ""
        outfile = ""

        next
    }

    print > outfile
}

' "${ANALYSIS_DIR}/representative-structures.pdb"

echo
echo "========================================"
echo " Clustering Analysis completed"
echo "========================================"
echo
echo "Setting:"
echo "  RMSD cutoff : ${CUTOFF} nm"
echo "  Frame skip  : ${SKIP}"
echo
echo "Output:"
echo "  ${ANALYSIS_DIR}/clustering-log.log"
echo "  ${ANALYSIS_DIR}/cluster-id.csv"
echo "  ${ANALYSIS_DIR}/cluster-index.ndx"
echo "  ${ANALYSIS_DIR}/cluster-size.xvg"
echo "  ${ANALYSIS_DIR}/representative-structures.pdb"
echo "  ${ANALYSIS_DIR}/cluster-○○.pdb"
echo " etc. "
echo
