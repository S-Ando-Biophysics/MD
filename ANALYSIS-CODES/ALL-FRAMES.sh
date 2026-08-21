#!/usr/bin/env bash

set -euo pipefail
export LC_ALL=C

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd -- "${SCRIPT_DIR}/.." && pwd)"

INPUT_DIR="${ROOT_DIR}/03_GROMACS/06_OUTPUT"

ANALYSIS_ROOT="${ROOT_DIR}/04_ANALYSIS"
OUTPUT_DIR="${ANALYSIS_ROOT}/00_ALL-FRAMES"

mkdir -p "$ANALYSIS_ROOT"
mkdir -p "$OUTPUT_DIR"

if ! command -v gmx >/dev/null 2>&1; then
    echo "ERROR: gmx command not found."
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
    echo "ERROR: Could not determine the number of frames from:"
    echo "  $XTC"
    echo
    echo "$CHECK_OUTPUT"
    exit 1
fi

if [ "$NUMBER_OF_FRAMES" -lt 1 ]; then
    echo "ERROR: Invalid number of frames:"
    echo "  $NUMBER_OF_FRAMES"
    exit 1
fi

LAST_FRAME_INDEX=$((NUMBER_OF_FRAMES - 1))

NZERO=${#LAST_FRAME_INDEX}

rm -f "${OUTPUT_DIR}"/frame*.pdb

printf "0\n" | \
gmx trjconv \
    -s "$PDB" \
    -f "$XTC" \
    -o "${OUTPUT_DIR}/frame.pdb" \
    -sep \
    -nzero "$NZERO"

NUMBER_OF_PDB_FILES=$(
    find "$OUTPUT_DIR" \
        -maxdepth 1 \
        -type f \
        -name "frame*.pdb" \
    | wc -l
)

