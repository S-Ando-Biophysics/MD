#!/usr/bin/env bash

set -euo pipefail
export LC_ALL=C

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd -- "${SCRIPT_DIR}/.." && pwd)"

INPUT_DIR="${ROOT_DIR}/03_GROMACS/06_OUTPUT"
ANALYSIS_ROOT="${ROOT_DIR}/04_ANALYSIS"

CPPTRAJ_INPUT="${SCRIPT_DIR}/H-BOND.in"

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

if ! command -v cpptraj >/dev/null 2>&1; then
    echo "ERROR: cpptraj command not found."
    exit 1
fi

if [ ! -f "$CPPTRAJ_INPUT" ]; then
    echo "ERROR: CPPTRAJ input file not found:"
    echo "  $CPPTRAJ_INPUT"
    exit 1
fi

if [ ! -d "$INPUT_DIR" ]; then
    echo "ERROR: Input directory not found:"
    echo "  $INPUT_DIR"
    exit 1
fi

XTC=$(
    find "$INPUT_DIR" \
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
    find "$INPUT_DIR" \
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
    echo "  $INPUT_DIR"
    exit 1
fi

if [ -z "${PDB:-}" ]; then
    echo "ERROR: PDB file not found in:"
    echo "  $INPUT_DIR"
    exit 1
fi

echo
echo "========================================"
echo " Hydrogen bond analysis "
echo "========================================"
echo

(
    cd "$ANALYSIS_DIR"

    cpptraj \
        -p "$PDB" \
        -y "$XTC" \
        -i "$CPPTRAJ_INPUT" \
        > cpptraj.log 2>&1
)

echo
echo "========================================"
echo " Hydrogen bond analysis completed "
echo "========================================"
echo
