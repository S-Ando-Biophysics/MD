#!/usr/bin/env bash

set -euo pipefail
export LC_ALL=C

PROC_DIR="04_PROCESSING-GROMACS"

TPR="03_GROMACS/04_MD/md.tpr"
XTC="03_GROMACS/04_MD/md.xtc"
INDEX="${PROC_DIR}/index.ndx"

# ------------------------------------------------------------
# Check input files
# ------------------------------------------------------------

if [ ! -f "$TPR" ]; then
    echo "ERROR: $TPR not found."
    exit 1
fi

if [ ! -f "$XTC" ]; then
    echo "ERROR: $XTC not found."
    exit 1
fi

if [ ! -f "$INDEX" ]; then
    echo "ERROR: $INDEX not found."
    exit 1
fi

mkdir -p "$PROC_DIR"

# ------------------------------------------------------------
# Select group
# ------------------------------------------------------------

echo
echo "========================================"
echo " Groups in $INDEX"
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

    read -rp "Select target group number: " TARGET_GROUP

    if ! [[ "$TARGET_GROUP" =~ ^[0-9]+$ ]]; then
        echo "ERROR: Please enter a group number."
        continue
    fi

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
echo "Selected group:"
echo "  Number : $TARGET_GROUP"
echo "  Name   : $TARGET_NAME"
echo

# ------------------------------------------------------------
# Step 1: PBC whole
# ------------------------------------------------------------

echo "========================================"
echo " Step 1/4 : PBC whole "
echo "========================================"

printf "0\n" | \
gmx trjconv \
    -s "$TPR" \
    -f "$XTC" \
    -n "$INDEX" \
    -o "${PROC_DIR}/md_whole.xtc" \
    -pbc whole

# ------------------------------------------------------------
# Step 2: PBC cluster
# ------------------------------------------------------------

echo
echo "========================================"
echo " Step 2/4 : PBC cluster "
echo "========================================"

printf "%s\n0\n" "$TARGET_GROUP" | \
gmx trjconv \
    -s "$TPR" \
    -f "${PROC_DIR}/md_whole.xtc" \
    -n "$INDEX" \
    -o "${PROC_DIR}/md_cluster.xtc" \
    -pbc cluster

# ------------------------------------------------------------
# Step 3: Center
# ------------------------------------------------------------

echo
echo "========================================"
echo " Step 3/4 : Center "
echo "========================================"

printf "%s\n0\n" "$TARGET_GROUP" | \
gmx trjconv \
    -s "$TPR" \
    -f "${PROC_DIR}/md_cluster.xtc" \
    -n "$INDEX" \
    -o "${PROC_DIR}/md_center.xtc" \
    -center \
    -pbc mol \
    -ur compact

# ------------------------------------------------------------
# Step 4: Fit rotation + translation
# ------------------------------------------------------------

echo
echo "========================================"
echo " Step 4/4 : Fit "
echo "========================================"

printf "%s\n%s\n0\n" \
    "$TARGET_GROUP" \
    "$TARGET_GROUP" | \
gmx trjconv \
    -s "$TPR" \
    -f "${PROC_DIR}/md_center.xtc" \
    -n "$INDEX" \
    -o "${PROC_DIR}/md_fit.xtc" \
    -center \
    -fit rot+trans

# ------------------------------------------------------------
# Finish
# ------------------------------------------------------------

echo
echo "========================================"
echo " Finished "
echo "========================================"
echo

echo "Target group : $TARGET_GROUP ($TARGET_NAME)"
echo

echo "Output files:"
echo "  ${PROC_DIR}/md_whole.xtc"
echo "  ${PROC_DIR}/md_cluster.xtc"
echo "  ${PROC_DIR}/md_center.xtc"
echo "  ${PROC_DIR}/md_fit.xtc"
echo