#!/usr/bin/env bash

set -euo pipefail
export LC_ALL=C

source ~/miniconda3/etc/profile.d/conda.sh

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd -- "${SCRIPT_DIR}/.." && pwd)"

ANALYSIS_ROOT="${ROOT_DIR}/04_ANALYSIS"
INPUT_DIR="${ANALYSIS_ROOT}/00_ALL-FRAMES"

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

OUTPUT_DIR="${ANALYSIS_ROOT}/${ANALYSIS_NUMBER}_FPOCKET-R"

mkdir -p "$OUTPUT_DIR"

WORK_DIR="$(mktemp -d)"

cleanup() {
    rm -rf "$WORK_DIR"
}

trap cleanup EXIT

conda activate fpocketR

echo
echo "========================================"
echo " Analysis using fpocketR "
echo "========================================"
echo
echo "Input:"
echo "  $INPUT_DIR"
echo
echo "Output:"
echo "  $OUTPUT_DIR"
echo

CURRENT_FILE=0

while IFS= read -r -d '' PDB_FILE; do

    CURRENT_FILE=$((CURRENT_FILE + 1))

    PDB_NAME="$(basename "$PDB_FILE")"
    PDB_STEM="${PDB_NAME%.pdb}"

    echo
    echo "------------------------------------------------------------"
    echo " ${CURRENT_FILE}/${NUMBER_OF_PDB_FILES} : ${PDB_NAME}"
    echo "------------------------------------------------------------"

    rm -rf "${WORK_DIR:?}/"*

    cp "$PDB_FILE" "${WORK_DIR}/${PDB_NAME}"

    mkdir -p "${WORK_DIR}/result"

    (
        cd "$WORK_DIR"

        python -m fpocketR \
            -pdb "$PDB_NAME" \
            -o result \
            -y
    )

    RESULT_DIR="${WORK_DIR}/result/${PDB_STEM}_clean_out"

    if [ ! -d "$RESULT_DIR" ]; then
        echo
        echo "ERROR: fpocketR output directory not found:"
        echo "  $RESULT_DIR"
        exit 1
    fi

    find "$RESULT_DIR" \
        -type f \
        ! -name "*info.txt" \
        ! -name "*pocket_characteristics.csv" \
        ! -name "*real_sphere.pdb" \
        ! -name "*png" \
        -delete

    find "$RESULT_DIR" \
        -mindepth 1 \
        -depth \
        -type d \
        -empty \
        -delete

    mv "$RESULT_DIR" "$OUTPUT_DIR/"

done < <(
    find "$INPUT_DIR" \
        -maxdepth 1 \
        -type f \
        -name "*.pdb" \
        -print0 \
    | sort -z
)

conda deactivate

echo
echo "========================================"
echo " fpocketR analysis completed "
echo "========================================"
echo
echo "Processed PDB files: $NUMBER_OF_PDB_FILES"
echo
echo "Output:"
echo "  $OUTPUT_DIR"
echo