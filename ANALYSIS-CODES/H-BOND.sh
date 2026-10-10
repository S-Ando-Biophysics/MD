#!/usr/bin/env bash

set -euo pipefail
export LC_ALL=C

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd -- "${SCRIPT_DIR}/.." && pwd)"

# ------------------------------------------------------------
# Analysis settings
# ------------------------------------------------------------

ANALYSIS_ROOT="${ROOT_DIR}/05_ANALYSIS"
CPPTRAJ_INPUT="${SCRIPT_DIR}/H-BOND.in"

# ------------------------------------------------------------
# Select post-processing method
# ------------------------------------------------------------

echo
echo "========================================"
echo " Hydrogen bond analysis "
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
# Check cpptraj
# ------------------------------------------------------------

if ! command -v cpptraj >/dev/null 2>&1; then
    echo
    echo "ERROR: cpptraj command not found."
    exit 1
fi

# ------------------------------------------------------------
# Check input files
# ------------------------------------------------------------

if [ ! -f "$CPPTRAJ_INPUT" ]; then
    echo
    echo "ERROR: CPPTRAJ input file not found:"
    echo "  $CPPTRAJ_INPUT"
    exit 1
fi

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

ANALYSIS_DIR="${ANALYSIS_ROOT}/${ANALYSIS_NUMBER}_HBOND"

mkdir -p "$ANALYSIS_DIR"

# ------------------------------------------------------------
# Hydrogen bond analysis
# ------------------------------------------------------------

echo
echo "========================================"
echo " Hydrogen bond analysis "
echo "========================================"
echo

echo "Post-processing method:"
echo "  $METHOD_NAME"
echo

echo "Input:"
echo "  PDB        : $PDB"
echo "  XTC        : $XTC"
echo "  CPPTRAJ IN : $CPPTRAJ_INPUT"
echo

echo "Output:"
echo "  $ANALYSIS_DIR"
echo

(
    cd "$ANALYSIS_DIR"

    cpptraj \
        -p "$PDB" \
        -y "$XTC" \
        -i "$CPPTRAJ_INPUT" \
        > cpptraj.log 2>&1
)

# ------------------------------------------------------------
# Finish
# ------------------------------------------------------------

echo
echo "========================================"
echo " Hydrogen bond analysis completed "
echo "========================================"
echo

echo "Post-processing method:"
echo "  $METHOD_NAME"
echo

echo "Output:"
echo "  ${ANALYSIS_DIR}"
echo "  ${ANALYSIS_DIR}/cpptraj.log"
echo