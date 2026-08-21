#!/usr/bin/env bash

set -euo pipefail
export LC_ALL=C

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd -- "${SCRIPT_DIR}/.." && pwd)"

ANALYSIS_ROOT="${ROOT_DIR}/04_ANALYSIS"
INPUT_DIR="${ANALYSIS_ROOT}/00_ALL-FRAMES"

if ! command -v find_pair >/dev/null 2>&1; then
    echo "ERROR: find_pair command not found."
    exit 1
fi

if ! command -v analyze >/dev/null 2>&1; then
    echo "ERROR: analyze command not found."
    exit 1
fi

if [ ! -d "$INPUT_DIR" ]; then
    echo "ERROR: Input directory not found:"
    echo "  $INPUT_DIR"
    exit 1
fi

NUMBER_OF_PDB_FILES=$(
    find "$INPUT_DIR" \
        -maxdepth 1 \
        -type f \
        -name "*.pdb" \
    | wc -l
)

if [ "$NUMBER_OF_PDB_FILES" -eq 0 ]; then
    echo "ERROR: No PDB files found in:"
    echo "  $INPUT_DIR"
    exit 1
fi

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

OUTPUT_DIR="${ANALYSIS_ROOT}/${ANALYSIS_NUMBER}_3DNA"

mkdir -p "$OUTPUT_DIR"

WORK_DIR="${OUTPUT_DIR}/tmp"

rm -rf "$WORK_DIR"
mkdir -p "$WORK_DIR"

echo
echo "========================================"
echo " Analysis using 3DNA "
echo "========================================"
echo

CURRENT_FILE=0

while IFS= read -r -d '' PDB_FILE; do

    CURRENT_FILE=$((CURRENT_FILE + 1))

    PDB_NAME="$(basename "$PDB_FILE")"
    PDB_STEM="${PDB_NAME%.pdb}"

    echo
    echo "----------------------------------------"
    echo " ${CURRENT_FILE}/${NUMBER_OF_PDB_FILES} : ${PDB_NAME}"
    echo "----------------------------------------"

    rm -rf "${WORK_DIR:?}/"*

    cp "$PDB_FILE" "${WORK_DIR}/${PDB_NAME}"

    (
        cd "$WORK_DIR"

        find_pair "$PDB_NAME" | analyze
    )

    OUT_FILE="${WORK_DIR}/${PDB_STEM}.out"

    if [ ! -f "$OUT_FILE" ]; then
        echo "ERROR: 3DNA output file was not generated:"
        echo "  ${PDB_STEM}.out"
        exit 1
    fi

    mv "$OUT_FILE" "$OUTPUT_DIR/"

done < <(
    find "$INPUT_DIR" \
        -maxdepth 1 \
        -type f \
        -name "*.pdb" \
        -print0 \
    | sort -z
)

rm -rf "$WORK_DIR"

find "$OUTPUT_DIR" \
    -maxdepth 1 \
    -type f \
    ! -name "*.out" \
    -delete

NUMBER_OF_OUT_FILES=$(
    find "$OUTPUT_DIR" \
        -maxdepth 1 \
        -type f \
        -name "*.out" \
    | wc -l
)

echo
echo "========================================"
echo " 3DNA analysis completed "
echo "========================================"
echo
echo "Input PDB files: $NUMBER_OF_PDB_FILES"
echo
echo "Output files: $NUMBER_OF_OUT_FILES"
echo
echo "Output:"
echo "  $OUTPUT_DIR"
echo