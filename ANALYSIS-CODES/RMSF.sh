#!/usr/bin/env bash

set -euo pipefail
export LC_ALL=C

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd -- "${SCRIPT_DIR}/.." && pwd)"

OUTPUT_DIR="${ROOT_DIR}/03_GROMACS/06_OUTPUT"

ANALYSIS_ROOT="${ROOT_DIR}/04_ANALYSIS"

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

ANALYSIS_DIR="${ANALYSIS_ROOT}/${ANALYSIS_NUMBER}_RMSF"

mkdir -p "$ANALYSIS_DIR"

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
echo " RMSF analysis"
echo "========================================"
echo

printf "0\n" | \
gmx rmsf \
    -s "$PDB" \
    -f "$XTC" \
    -res \
    -o "${ANALYSIS_DIR}/rmsf.xvg"

awk '
BEGIN {
    print "Residue,RMSF [nm],RMSF [A]"
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
' "${ANALYSIS_DIR}/rmsf.xvg" \
> "${ANALYSIS_DIR}/rmsf.csv"

echo
echo "========================================"
echo " RMSF analysis completed"
echo "========================================"
echo
echo "Output:"
echo "  ${ANALYSIS_DIR}/rmsf.xvg"
echo "  ${ANALYSIS_DIR}/rmsf.csv"
echo
