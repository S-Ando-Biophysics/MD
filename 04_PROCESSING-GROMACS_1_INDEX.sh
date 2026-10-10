#!/usr/bin/env bash

set -euo pipefail
export LC_ALL=C

PROC_DIR="04_PROCESSING-GROMACS"
TPR="03_GROMACS/04_MD/md.tpr"
INDEX="${PROC_DIR}/index.ndx"

# ------------------------------------------------------------
# Check input file
# ------------------------------------------------------------

if [ ! -f "$TPR" ]; then
    echo "ERROR: TPR file not found:"
    echo "  $TPR"
    exit 1
fi

# ------------------------------------------------------------
# Create processing directory
# ------------------------------------------------------------

mkdir -p "$PROC_DIR"

# ------------------------------------------------------------
# Create index
# ------------------------------------------------------------

echo
echo "========================================"
echo " Create index "
echo "========================================"
echo

echo "Input:"
echo "  TPR   : $TPR"
echo

echo "Output:"
echo "  INDEX : $INDEX"
echo

gmx make_ndx \
    -f "$TPR" \
    -o "$INDEX"

# l
# ri 1-○○
# name ○○ Nuc
# q

echo
echo "========================================"
echo " Index creation completed "
echo "========================================"
echo

echo "Created:"
echo "  $INDEX"
echo