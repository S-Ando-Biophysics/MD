#!/usr/bin/env bash

set -euo pipefail
export LC_ALL=C


SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd -- "${SCRIPT_DIR}/.." && pwd)"


# ------------------------------------------------------------
# Analysis directory
# ------------------------------------------------------------

ANALYSIS_ROOT="${ROOT_DIR}/05_ANALYSIS"


# ------------------------------------------------------------
# Select post-processing method
# ------------------------------------------------------------

echo
echo "========================================"
echo " Clustering Analysis "
echo "========================================"
echo

echo "Select trajectory post-processing method:"
echo
echo "  1) GROMACS"
echo "  2) AmberTools"
echo

read -rp "Selection [1/2]: " PROC_METHOD


case "$PROC_METHOD" in

    1)
        METHOD_NAME="GROMACS"

        PDB="${ROOT_DIR}/04_PROCESSING-GROMACS/md_target.pdb"
        XTC="${ROOT_DIR}/04_PROCESSING-GROMACS/md_fit_target.xtc"
        INDEX="${ROOT_DIR}/04_PROCESSING-GROMACS/index.ndx"
        ;;

    2)
        METHOD_NAME="AmberTools"

        PDB="${ROOT_DIR}/04_PROCESSING-AMBER/md_target.pdb"
        XTC="${ROOT_DIR}/04_PROCESSING-AMBER/md_fit_target.xtc"
        INDEX="${ROOT_DIR}/04_PROCESSING-AMBER/index.ndx"
        ;;

    *)
        echo
        echo "ERROR: Invalid selection."
        echo "Please select 1 or 2."
        exit 1
        ;;

esac


# ------------------------------------------------------------
# Check GROMACS
# ------------------------------------------------------------

if ! command -v gmx >/dev/null 2>&1; then
    echo
    echo "ERROR: gmx command not found."
    exit 1
fi


# ------------------------------------------------------------
# Check input files
# ------------------------------------------------------------

if [ ! -f "$PDB" ]; then
    echo
    echo "ERROR: PDB file not found:"
    echo "  $PDB"
    exit 1
fi


if [ ! -f "$XTC" ]; then
    echo
    echo "ERROR: XTC file not found:"
    echo "  $XTC"
    exit 1
fi


if [ ! -f "$INDEX" ]; then
    echo
    echo "ERROR: index file not found:"
    echo "  $INDEX"
    exit 1
fi


# ------------------------------------------------------------
# Create analysis directory
# ------------------------------------------------------------

mkdir -p "$ANALYSIS_ROOT"


MAX_NUMBER=$(
    find "$ANALYSIS_ROOT" \
        -mindepth 1 \
        -maxdepth 1 \
        -type d \
        -printf "%f\n" \
    | awk '
        /^[0-9]+_/ {
            split($0, parts, "_")
            number = parts[1] + 0

            if (number > max) {
                max = number
            }
        }

        END {
            print max + 0
        }
    '
)


NEXT_NUMBER=$((MAX_NUMBER + 1))


printf -v ANALYSIS_NUMBER "%02d" "$NEXT_NUMBER"


ANALYSIS_DIR="${ANALYSIS_ROOT}/${ANALYSIS_NUMBER}_CLUSTER-GROMOS"


mkdir -p "$ANALYSIS_DIR"


# ------------------------------------------------------------
# RMSD cutoff
# ------------------------------------------------------------

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


# ------------------------------------------------------------
# Frame skip
# ------------------------------------------------------------

while true; do

    read -rp "Analyze every N-th frame: " SKIP

    if [[ "$SKIP" =~ ^[1-9][0-9]*$ ]]; then
        break
    fi

    echo "ERROR: Please enter a positive integer."
    echo "Example: 10"

done


# ------------------------------------------------------------
# Clustering
# ------------------------------------------------------------

echo
echo "========================================"
echo " Clustering Analysis "
echo "========================================"
echo

echo "Post-processing method:"
echo "  $METHOD_NAME"
echo

echo "Input:"
echo "  PDB   : $PDB"
echo "  XTC   : $XTC"
echo "  INDEX : $INDEX"
echo

echo "Setting:"
echo "  Method      : GROMOS"
echo "  RMSD cutoff : ${CUTOFF} nm"
echo "  Frame skip  : ${SKIP}"
echo

echo "Output:"
echo "  $ANALYSIS_DIR"
echo


gmx cluster \
    -s "$PDB" \
    -f "$XTC" \
    -n "$INDEX" \
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


# ------------------------------------------------------------
# Convert cluster ID XVG to CSV
# ------------------------------------------------------------

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


# ------------------------------------------------------------
# Split representative structures
# ------------------------------------------------------------

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


# ------------------------------------------------------------
# Finish
# ------------------------------------------------------------

echo
echo "========================================"
echo " Clustering Analysis completed "
echo "========================================"
echo

echo "Post-processing method:"
echo "  $METHOD_NAME"
echo

echo "Input:"
echo "  PDB   : $PDB"
echo "  XTC   : $XTC"
echo "  INDEX : $INDEX"
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
echo "  etc."
echo