#!/bin/bash

set -euo pipefail

cd 03_GROMACS/05_PROC-GROMACS

TPR="../04_MD/md.tpr"
XTC="../04_MD/md.xtc"
INDEX="index.ndx"

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


# ------------------------------------------------------------
# Select group
# ------------------------------------------------------------

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
#
# Selection:
#   1st = System (0)
# ------------------------------------------------------------

echo "========================================"
echo " Step 1/4 : PBC whole "
echo "========================================"

printf "0\n" | \
gmx trjconv \
    -s "$TPR" \
    -f "$XTC" \
    -n "$INDEX" \
    -o md_whole.xtc \
    -pbc whole


# ------------------------------------------------------------
# Step 2: PBC cluster
#
# Selection:
#   1st = selected target group
#   2nd = System (0)
# ------------------------------------------------------------

echo
echo "========================================"
echo " Step 2/4 : PBC cluster "
echo "========================================"

printf "%s\n0\n" "$TARGET_GROUP" | \
gmx trjconv \
    -s "$TPR" \
    -f md_whole.xtc \
    -n "$INDEX" \
    -o md_cluster.xtc \
    -pbc cluster


# ------------------------------------------------------------
# Step 3. Center
#
# Selection:
#   1st = selected target group
#   2nd = System (0)
# ------------------------------------------------------------

echo
echo "========================================"
echo " Step 3/4 : Center "
echo "========================================"

printf "%s\n0\n" "$TARGET_GROUP" | \
gmx trjconv \
    -s "$TPR" \
    -f md_cluster.xtc \
    -n "$INDEX" \
    -o md_center.xtc \
    -center \
    -pbc mol \
    -ur compact


# ------------------------------------------------------------
# Step 4. Fit rotation + translation
#
# Selection:
#   1st = selected target group
#   2nd = selected target group
#   3rd = System (0)
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
    -f md_center.xtc \
    -n "$INDEX" \
    -o md_fit.xtc \
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
echo "  md_whole.xtc"
echo "  md_cluster.xtc"
echo "  md_center.xtc"
echo "  md_fit.xtc"
echo

cd ../..
