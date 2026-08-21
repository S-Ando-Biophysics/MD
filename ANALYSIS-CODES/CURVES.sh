#!/usr/bin/env bash

set -euo pipefail
export LC_ALL=C


# ============================================================
# User settings
# ============================================================

# Directory where Curves+ is installed
CURVES_DIR="/path/to/Curves+"

# Residue numbers of 1st strand (5'→3')
STRAND_1="1:12"

# Residue numbers of 2nd strand (3'→5')
STRAND_2="24:13"


# ============================================================


SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd -- "${SCRIPT_DIR}/.." && pwd)"

ANALYSIS_ROOT="${ROOT_DIR}/04_ANALYSIS"
INPUT_DIR="${ANALYSIS_ROOT}/00_ALL-FRAMES"

if [ ! -d "$CURVES_DIR" ]; then
    echo "ERROR: Curves+ directory not found:"
    echo "  $CURVES_DIR"
    exit 1
fi

if [ ! -x "${CURVES_DIR}/Cur+" ]; then
    echo "ERROR: Cur+ executable not found:"
    echo "  ${CURVES_DIR}/Cur+"
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

OUTPUT_DIR="${ANALYSIS_ROOT}/${ANALYSIS_NUMBER}_CURVES"

mkdir -p "$OUTPUT_DIR"

echo
echo "========================================"
echo " Analysis using Curves+ "
echo "========================================"
echo

CURRENT_FILE=0

while IFS= read -r -d '' PDB_FILE; do

    CURRENT_FILE=$((CURRENT_FILE + 1))

    PDB_NAME="$(basename "$PDB_FILE")"
    PDB_STEM="${PDB_NAME%.pdb}"

    TEMP_STEM="CURVES_TMP_${$}_${CURRENT_FILE}"
    TEMP_PDB="${TEMP_STEM}.pdb"
    TEMP_OUTPUT="${TEMP_STEM}-Output"

    echo
    echo "----------------------------------------"
    echo " ${CURRENT_FILE}/${NUMBER_OF_PDB_FILES} : ${PDB_NAME}"
    echo "----------------------------------------"

    cp "$PDB_FILE" "${CURVES_DIR}/${TEMP_PDB}"

    if ! (
        cd "$CURVES_DIR"

        ./Cur+ << EOF
&inp file=${TEMP_PDB},lis=${TEMP_OUTPUT},lib=standard, &end
2 1 -1 0 0
${STRAND_1}
${STRAND_2}
EOF
    ); then

        rm -f \
            "${CURVES_DIR}/${TEMP_PDB}" \
            "${CURVES_DIR}/${TEMP_OUTPUT}.lis" \
            "${CURVES_DIR}/${TEMP_OUTPUT}.cda" \
            "${CURVES_DIR}/${TEMP_OUTPUT}_B.pdb" \
            "${CURVES_DIR}/${TEMP_OUTPUT}_X.pdb" \
            "${CURVES_DIR}/${TEMP_OUTPUT}_R.pdb"

        echo "ERROR: Curves+ analysis failed:"
        echo "  $PDB_NAME"
        exit 1
    fi

    LIS_FILE="${CURVES_DIR}/${TEMP_OUTPUT}.lis"

    if [ ! -f "$LIS_FILE" ]; then
        rm -f \
            "${CURVES_DIR}/${TEMP_PDB}" \
            "${CURVES_DIR}/${TEMP_OUTPUT}.cda" \
            "${CURVES_DIR}/${TEMP_OUTPUT}_B.pdb" \
            "${CURVES_DIR}/${TEMP_OUTPUT}_X.pdb" \
            "${CURVES_DIR}/${TEMP_OUTPUT}_R.pdb"

        echo "ERROR: Curves+ output file was not generated:"
        echo "  ${PDB_STEM}-Output.lis"
        exit 1
    fi

    mv \
        "$LIS_FILE" \
        "${OUTPUT_DIR}/${PDB_STEM}-Output.lis"

    rm -f \
        "${CURVES_DIR}/${TEMP_PDB}" \
        "${CURVES_DIR}/${TEMP_OUTPUT}.cda" \
        "${CURVES_DIR}/${TEMP_OUTPUT}_B.pdb" \
        "${CURVES_DIR}/${TEMP_OUTPUT}_X.pdb" \
        "${CURVES_DIR}/${TEMP_OUTPUT}_R.pdb"

done < <(
    find "$INPUT_DIR" \
        -maxdepth 1 \
        -type f \
        -name "*.pdb" \
        -print0 \
    | sort -z
)

NUMBER_OF_LIS_FILES=$(
    find "$OUTPUT_DIR" \
        -maxdepth 1 \
        -type f \
        -name "*-Output.lis" \
    | wc -l
)

echo
echo "========================================"
echo " Curves+ analysis completed "
echo "========================================"
echo
echo "Input PDB files: $NUMBER_OF_PDB_FILES"
echo
echo "Output files: $NUMBER_OF_LIS_FILES"
echo
echo "Output:"
echo "  $OUTPUT_DIR"
echo
