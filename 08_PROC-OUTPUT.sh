#!/usr/bin/env bash

set -euo pipefail

PROC_DIR="03_GROMACS/05_PROC-GROMACS"
OUTPUT_DIR="03_GROMACS/06_OUTPUT"

TPR="03_GROMACS/04_MD/md.tpr"
INDEX="${PROC_DIR}/index.ndx"
XTC="${PROC_DIR}/md_fit.xtc"

if [ ! -f "$TPR" ]; then
    echo "ERROR: $TPR not found."
    exit 1
fi

if [ ! -f "$INDEX" ]; then
    echo "ERROR: $INDEX not found."
    exit 1
fi

if [ ! -f "$XTC" ]; then
    echo "ERROR: $XTC not found."
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

echo
echo "========================================"
echo " Groups in index.ndx"
echo "========================================"

awk '
/^[[:space:]]*\[/ {
    name = $0

    gsub(/^[[:space:]]*\[[[:space:]]*/, "", name)
    gsub(/[[:space:]]*\][[:space:]]*$/, "", name)

    printf "%3d : %s\n", n, name

    n++
}
' "$INDEX"

echo "========================================"
echo

N_GROUPS=$(awk '
/^[[:space:]]*\[/ {
    n++
}
END {
    print n
}
' "$INDEX")

while true; do

    read -rp "Select output group number: " TARGET_GROUP

    # Integer check
    if ! [[ "$TARGET_GROUP" =~ ^[0-9]+$ ]]; then
        echo "ERROR: Please enter a valid group number."
        continue
    fi

    # Range check
    if [ "$TARGET_GROUP" -ge "$N_GROUPS" ]; then
        echo "ERROR: Group $TARGET_GROUP does not exist."
        continue
    fi

    break

done

TARGET_NAME=$(awk -v target="$TARGET_GROUP" '
/^[[:space:]]*\[/ {

    name = $0

    gsub(/^[[:space:]]*\[[[:space:]]*/, "", name)
    gsub(/[[:space:]]*\][[:space:]]*$/, "", name)

    if (n == target) {
        print name
        exit
    }

    n++
}
' "$INDEX")

echo
echo "Selected output group:"
echo "  Number : $TARGET_GROUP"
echo "  Name   : $TARGET_NAME"
echo

printf "%s\n" "$TARGET_GROUP" | \
gmx trjconv \
    -s "$TPR" \
    -f "$XTC" \
    -n "$INDEX" \
    -o "${OUTPUT_DIR}/md_fit_target.xtc"

printf "%s\n" "$TARGET_GROUP" | \
gmx trjconv \
    -s "$TPR" \
    -f "$XTC" \
    -n "$INDEX" \
    -o "${OUTPUT_DIR}/md_target.pdb" \
    -dump 0

echo
echo "========================================"
echo " Finished"
echo "========================================"
echo
echo "Selected group:"
echo "  $TARGET_GROUP : $TARGET_NAME"
echo
echo "Output files:"
echo "  ${OUTPUT_DIR}/md_fit_target.xtc"
echo "  ${OUTPUT_DIR}/md_target.pdb"
echo
