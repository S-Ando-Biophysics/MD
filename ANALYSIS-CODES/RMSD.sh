#!/usr/bin/env bash

set -euo pipefail
export LC_ALL=C

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd -- "${SCRIPT_DIR}/.." && pwd)"

# ------------------------------------------------------------
# Input files
# ------------------------------------------------------------

TPR="${ROOT_DIR}/03_GROMACS/04_MD/md.tpr"
XTC="${ROOT_DIR}/03_GROMACS/05_PROC-GROMACS/md_fit.xtc"
INDEX="${ROOT_DIR}/03_GROMACS/05_PROC-GROMACS/index.ndx"

ANALYSIS_ROOT="${ROOT_DIR}/04_ANALYSIS"

# ------------------------------------------------------------
# Check input files
# ------------------------------------------------------------

if [ ! -f "$TPR" ]; then
    echo "ERROR: TPR file not found:"
    echo "  $TPR"
    exit 1
fi

if [ ! -f "$XTC" ]; then
    echo "ERROR: XTC file not found:"
    echo "  $XTC"
    exit 1
fi

if [ ! -f "$INDEX" ]; then
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

ANALYSIS_DIR="${ANALYSIS_ROOT}/${ANALYSIS_NUMBER}_RMSD"

mkdir -p "$ANALYSIS_DIR"

# ------------------------------------------------------------
# RMSD analysis
# ------------------------------------------------------------

echo
echo "========================================"
echo " RMSD analysis "
echo "========================================"
echo

echo "Input:"
echo "  TPR   : $TPR"
echo "  XTC   : $XTC"
echo "  INDEX : $INDEX"
echo

gmx rms \
    -s "$TPR" \
    -f "$XTC" \
    -n "$INDEX" \
    -tu ns \
    -o "${ANALYSIS_DIR}/rmsd.xvg"

# ------------------------------------------------------------
# Convert XVG to CSV
# ------------------------------------------------------------

awk '
BEGIN {
    print "Time [ns],RMSD [nm],RMSD [A]"
}

/^[[:space:]]*#/ {
    next
}

/^[[:space:]]*@/ {
    next
}

NF >= 2 {
    print $1 "," $2 "," ($2 * 10)
}
' "${ANALYSIS_DIR}/rmsd.xvg" \
> "${ANALYSIS_DIR}/rmsd.csv"

# ------------------------------------------------------------
# Finish
# ------------------------------------------------------------

echo
echo "========================================"
echo " RMSD analysis completed"
echo "========================================"
echo
echo "Output:"
echo "  ${ANALYSIS_DIR}/rmsd.xvg"
echo "  ${ANALYSIS_DIR}/rmsd.csv"
echo
