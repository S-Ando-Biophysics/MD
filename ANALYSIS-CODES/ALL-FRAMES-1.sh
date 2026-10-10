#!/usr/bin/env bash

set -euo pipefail
export LC_ALL=C

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd -- "${SCRIPT_DIR}/.." && pwd)"

# ------------------------------------------------------------
# Analysis directory
# ------------------------------------------------------------

ANALYSIS_ROOT="${ROOT_DIR}/05_ANALYSIS"
OUTPUT_DIR="${ANALYSIS_ROOT}/00_ALL-FRAMES"

# ------------------------------------------------------------
# Select post-processing method
# ------------------------------------------------------------

echo
echo "========================================"
echo " Extract all frames "
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
        ;;

    2)
        METHOD_NAME="AmberTools"

        PDB="${ROOT_DIR}/04_PROCESSING-AMBER/md_target.pdb"
        XTC="${ROOT_DIR}/04_PROCESSING-AMBER/md_fit_target.xtc"
        ;;

    *)
        echo
        echo "ERROR: Invalid selection."
        echo "Please select 1 or 2."
        exit 1
        ;;

esac

# ------------------------------------------------------------
# Check command
# ------------------------------------------------------------

if ! command -v gmx >/dev/null 2>&1; then
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

# ------------------------------------------------------------
# Create output directory
# ------------------------------------------------------------

mkdir -p "$ANALYSIS_ROOT"
mkdir -p "$OUTPUT_DIR"

# ------------------------------------------------------------
# Determine number of frames
# ------------------------------------------------------------

CHECK_OUTPUT=$(
    gmx check -f "$XTC" 2>&1
)

NUMBER_OF_FRAMES=$(
    printf "%s\n" "$CHECK_OUTPUT" \
    | awk '
        $1 == "Step" && $2 ~ /^[0-9]+$/ {
            print $2
            exit
        }
    '
)

if [ -z "${NUMBER_OF_FRAMES:-}" ]; then

    NUMBER_OF_FRAMES=$(
        printf "%s\n" "$CHECK_OUTPUT" \
        | awk '
            $1 == "Coords" && $2 ~ /^[0-9]+$/ {
                print $2
                exit
            }
        '
    )

fi

if [ -z "${NUMBER_OF_FRAMES:-}" ]; then
    echo
    echo "ERROR: Could not determine the number of frames from:"
    echo "  $XTC"
    echo
    echo "$CHECK_OUTPUT"
    exit 1
fi

if [ "$NUMBER_OF_FRAMES" -lt 1 ]; then
    echo
    echo "ERROR: Invalid number of frames:"
    echo "  $NUMBER_OF_FRAMES"
    exit 1
fi

LAST_FRAME_INDEX=$((NUMBER_OF_FRAMES - 1))

NZERO=${#LAST_FRAME_INDEX}

# ------------------------------------------------------------
# Extract all frames
# ------------------------------------------------------------

echo
echo "========================================"
echo " Extract all frames "
echo "========================================"
echo

echo "Post-processing method:"
echo "  $METHOD_NAME"
echo

echo "Input:"
echo "  PDB : $PDB"
echo "  XTC : $XTC"
echo

echo "Frames:"
echo "  $NUMBER_OF_FRAMES"
echo

echo "Output:"
echo "  $OUTPUT_DIR"
echo

rm -f "${OUTPUT_DIR}"/frame*.pdb

printf "0\n" | \
gmx trjconv \
    -s "$PDB" \
    -f "$XTC" \
    -o "${OUTPUT_DIR}/frame.pdb" \
    -sep \
    -nzero "$NZERO"

# ------------------------------------------------------------
# Count output files
# ------------------------------------------------------------

NUMBER_OF_PDB_FILES=$(
    find "$OUTPUT_DIR" \
        -maxdepth 1 \
        -type f \
        -name "frame*.pdb" \
    | wc -l
)

# ------------------------------------------------------------
# Finish
# ------------------------------------------------------------

echo
echo "========================================"
echo " All-frame extraction completed "
echo "========================================"
echo

echo "Post-processing method:"
echo "  $METHOD_NAME"
echo

echo "Number of frames:"
echo "  $NUMBER_OF_PDB_FILES"
echo

echo "Output:"
echo "  $OUTPUT_DIR"
echo